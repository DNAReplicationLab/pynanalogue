//! # Pynanalogue (Python functions for Nanalogue (= Nucleic Acid Analogue))
//!
//! In Nanalogue, we process and analyse data associated with DNA/RNA molecules,
//! their alignments to reference genomes, modification information on them,
//! and other miscellaneous information from BAM files. We expose some of
//! these functions in a python package for usage by the python bioinformatics
//! community.
//!
//! ## Design Philosophy: Many Parameters vs Parameter Objects
//!
//! The exported Python functions in this module accept numerous individual parameters
//! (typically 15-20) rather than using parameter objects or configuration structs.
//! This design intentionally prioritizes Python user experience over Rust code elegance.
//! Most parameters are optional with sensible defaults, allowing Python users to write
//! clear, self-documenting code like `read_info(bam_path, min_seq_len=1000, mapq_filter=20)`
//! with full IDE autocomplete support. While this creates some code duplication in the
//! Rust implementation, it provides a familiar, low-friction API that matches common
//! patterns in scientific Python libraries. For a solo-maintained side project, this
//! tradeoff values shipping a usable Python package over internal Rust maintainability.
use nanalogue_core::{
    BamPreFilt as _, BamRcRecords, CurrRead, Error, F32Bw0and1, GenomicRegion, InputBam,
    InputBamBuilder, InputMods, InputModsBuilder, InputWindowingBuilder, OptionalTag, OrdPair,
    PathOrURLOrStdin, SeqDisplayOptions, SimulationConfig, ThresholdState, analysis,
    curr_reads_to_dataframe, nanalogue_indexed_bam_reader, nanalogue_indexed_bam_reader_from_url,
    peek as rust_peek, read_info as rust_read_info, reads_table as rust_reads_table,
    simulate_mod_bam as rust_simulate_mod_bam, window_reads as rust_window_reads,
};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use pyo3_polars::PyDataFrame;
use rust_htslib::bam::FetchDefinition;
use std::collections::HashSet;
use std::num::NonZeroU32;
use std::str::FromStr as _;
use url::{ParseError, Url};

// Avoid musl's default allocator due to lackluster performance
// https://nickb.dev/blog/default-musl-allocator-considered-harmful-to-performance
#[cfg(target_env = "musl")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Converts a `nanalogue_core::Error` to a `PyException` through `Display`
macro_rules! py_exception {
    ($a:expr) => {
        pyo3::exceptions::PyException::new_err($a.to_string())
    };
}

/// Converts a `nanalogue_core::Error` to a `PyIOError` through `Display`
macro_rules! py_io_error {
    ($a:expr) => {
        pyo3::exceptions::PyIOError::new_err($a.to_string())
    };
}

/// Converts a `nanalogue_core::Error` to a `PyValueError` through `Display`
macro_rules! py_value_error {
    ($a:expr) => {
        pyo3::exceptions::PyValueError::new_err($a.to_string())
    };
}

/// Parse input options and convert them into our `InputBam`, `InputMods` structs.
///
/// # Input Validation Strategy
///
/// This function does minimal validation at the Python wrapper level.
/// Detailed input validation and error messages are provided by the underlying
/// `nanalogue` Rust library. This avoids code duplication and ensures consistency
/// between the Rust CLI and Python bindings. Invalid inputs (e.g., malformed
/// region strings, invalid filter names) will produce descriptive error messages
/// from the `nanalogue` library when the builder methods are called.
#[expect(
    clippy::let_underscore_untyped,
    reason = "occasionally we will leave _ untyped"
)]
#[expect(
    clippy::too_many_arguments,
    clippy::fn_params_excessive_bools,
    reason = "See module-level doc 'Design Philosophy: Many Parameters vs Parameter Objects'"
)]
fn parse_input_options(
    bam_path: &str,
    treat_as_url: bool,
    min_seq_len: u64,
    min_align_len: i64,
    read_ids: HashSet<String>,
    threads: u8,
    include_zero_len: bool,
    read_filter: &str,
    sample_fraction: f32,
    mapq_filter: u8,
    exclude_mapq_unavail: bool,
    region: &str,
    full_region: bool,
    tag: &str,
    mod_strand: &str,
    min_mod_qual: u8,
    reject_mod_qual_non_inclusive: (u8, u8),
    trim_read_ends_mod: usize,
    base_qual_filter_mod: u8,
    mod_region: &str,
) -> PyResult<(InputBam, InputMods<OptionalTag>)> {
    // Guard against zero-length reads until core is hardened
    if include_zero_len {
        return Err(py_value_error!(
            "include_zero_len=True is not yet supported due to potential crashes in the underlying library. \
             Please filter zero-length reads upstream or wait for this feature to be stabilized."
        ));
    }

    let mut bam_builder = InputBamBuilder::default();
    let _ = bam_builder
        .bam_path({
            if treat_as_url {
                PathOrURLOrStdin::URL(
                    Url::parse(bam_path).map_err(|e: ParseError| py_value_error!(e))?,
                )
            } else {
                PathOrURLOrStdin::Path(bam_path.into())
            }
        })
        .min_seq_len(min_seq_len)
        .threads(
            NonZeroU32::new(threads.into())
                .ok_or(py_value_error!("threads must be a positive integer"))?,
        )
        .include_zero_len(include_zero_len)
        .read_filter({
            // Trim whitespace from comma-separated tokens
            read_filter
                .split(',')
                .map(str::trim)
                .collect::<Vec<_>>()
                .join(",")
        })
        .sample_fraction(F32Bw0and1::new(sample_fraction).map_err(|e| py_value_error!(e))?)
        .mapq_filter(mapq_filter)
        .exclude_mapq_unavail(exclude_mapq_unavail)
        .region(region.into())
        .full_region(full_region);
    // Validate and set min_align_len based on value
    match min_align_len {
        ..0 => {
            return Err(py_value_error!(format!(
                "min_align_len must be non-negative, got {}",
                min_align_len
            )));
        }
        0 => {} // Don't set if zero (default behavior)
        v => {
            let _ = bam_builder.min_align_len(v);
        }
    }
    let _: Option<&mut _> = (!read_ids.is_empty()).then(|| bam_builder.read_id_set(read_ids));

    let bam = bam_builder.build().map_err(|e| py_value_error!(e))?;

    let mut mods_builder = InputModsBuilder::<OptionalTag>::default();
    let _ = mods_builder
        .mod_strand(mod_strand.into())
        .mod_prob_filter({
            let low = reject_mod_qual_non_inclusive.0;
            let high = reject_mod_qual_non_inclusive.1;
            match high.checked_sub(low) {
                None => {
                    return Err(py_value_error!(
                        "for rejecting mod quals, please set low < high"
                    ));
                }
                Some(0 | 1) => ThresholdState::GtEq(min_mod_qual),
                _ => {
                    #[expect(
                        clippy::arithmetic_side_effects,
                        reason = "we check low < high - 1 so no chance of over or underflow"
                    )]
                    let ord_pair = OrdPair::<u8>::try_from((low + 1, high - 1))
                        .map_err(|e| py_value_error!(e))?;
                    ThresholdState::Both((min_mod_qual, ord_pair))
                }
            }
        })
        .trim_read_ends_mod(trim_read_ends_mod)
        .base_qual_filter_mod(base_qual_filter_mod)
        .mod_region(mod_region.into());
    if !tag.is_empty() {
        let _ = mods_builder.tag(OptionalTag::from_str(tag).map_err(|e| py_value_error!(e))?);
    }
    let mods = mods_builder.build().map_err(|e| py_value_error!(e))?;

    Ok((bam, mods))
}

/// Load BAM data from a local file or a URL; fetch only the region if specified.
/// Needs an associated BAM index.
fn load_bam(bam: InputBam) -> PyResult<rust_htslib::bam::IndexedReader> {
    let reader = match (bam.region, bam.bam_path) {
        (Some(v), PathOrURLOrStdin::Path(w)) => nanalogue_indexed_bam_reader(
            &w,
            (&v).try_into().map_err(|e: Error| py_value_error!(e))?,
        ),
        (None, PathOrURLOrStdin::Path(w)) => nanalogue_indexed_bam_reader(&w, FetchDefinition::All),
        (Some(v), PathOrURLOrStdin::URL(w)) => nanalogue_indexed_bam_reader_from_url(
            &w,
            (&v).try_into().map_err(|e: Error| py_value_error!(e))?,
        ),
        (None, PathOrURLOrStdin::URL(w)) => {
            nanalogue_indexed_bam_reader_from_url(&w, FetchDefinition::All)
        }
        _ => {
            return Err(py_value_error!(
                "Stdin input not supported for indexed BAM reading"
            ));
        }
    }
    .map_err(|e: Error| py_io_error!(e))?;
    Ok(reader)
}

/// Produces bytes which can be decoded to JSON format,
/// with one JSON record per BAM record with information per record such as
/// alignment length, sequence length, read id, mod counts etc.
///
/// Sets various options through builder functions before running
/// the function and capturing the output as a stream of bytes
/// which can be decoded to JSON (see the Example output section below).
/// Runs `nanalogue_core::read_info::run`.
///
/// This function can be used to count how many reads are above a threshold
/// length or modification count or are primary mappings.
/// The function can also be used to analyze relationships such as alignment
/// length vs basecalled length. The arguments can be used to filter
/// the BAM data (e.g. passing through a specific region etc.).
///
/// # Args
/// bam_path (str): Path to the BAM file. Must be associated with a BAM index.
/// treat_as_url (optional, bool): If True, treat `bam_path` as a URL, default False.
/// min_seq_len (optional, int): Only retain sequences above this length, default 0.
/// min_align_len (optional, int): Only retain sequences with an alignment length above this
///     value. Defaults to unused.
/// read_ids (optional, set of str): Only retrieve these read ids, defaults to unused.
/// threads (optional, int): Number of threads used in some aspects of program execution, defaults to 2.
/// include_zero_len (optional, bool): Include sequences of zero length. WARNING: our program
///     may crash if you do this. Defaults to False. Helps to check if sequences of zero length
///     exist in our BAM file.
/// read_filter (optional, str): Comma-separated sequence of one to many of the following
///     strings: primary_forward, primary_reverse, secondary_forward, secondary_reverse,
///     supplementary_forward, supplementary_reverse, unmapped. If specified, only reads
///     with a mapping belonging to this set are retained. Defaults to no filter.
/// sample_fraction (optional, float): Set to between 0 and 1 to subsample BAM file.
///     WARNING: seeds are not set, so you may get a new set of reads every time.
///     WARNING: we sample every read with the given probability, so the total number
///         of reads fluctuates according to standard counting statistics.
/// mapq_filter (optional, int): Exclude reads with mapping quality below this number.
///     defaults to unused.
/// exclude_mapq_unavail (optional, bool): Exclude reads where mapping quality is unavailable.
///     defaults to false.
/// region (optional, str): Only include reads with at least one mapped base from this region.
///     Use the format "contig", "contig:start-", or "contig:start-end". These are 0-based,
///     half-open intervals. Defaults to read entire BAM file. Can be used in combination
///     with `mod_region`.
/// full_region (optional, bool): Only include reads if they pass through the region above
///     in full. Defaults to false.
/// tag (optional, str): If set, only this type of modification is processed. Input is a
///     string type, for example a single-letter code "m", a number as a string "76792" etc.
///     Defaults to processing all modifications.
/// mod_strand (optional, str): Set this to `bc` or `bc_comp` to retrieve information
///     about mods only from the basecalled strand or only from its complement.
///     Some sequencing technologies like `PacBio` or `ONT duplex` record mod information
///     both on a strand and its complement. It may be useful in some scenarios to
///     separate this information. Defaults to not filter.
/// min_mod_qual (optional, int): Set to a number 0-255. Reject modification
///     calls whose probability is below this value (0, 255 correspond to a
///     probability of 0 and 1 respectively). Defaults to 0.
/// reject_mod_qual_non_inclusive (optional, (int, int)): Reject modification
///     calls whose probability is such that int_low < prob < int_high.
///     Set both to a number between 0-255 and such that the first entry is <=
///     the second (if they are equal, no filtering is performed). Defaults
///     to no filtering. Also see comments under `min_mod_qual`.
/// trim_read_ends_mod (optional, int): Reject modification information
///     within so many bp of either end of the read. Defaults to 0.
/// base_qual_filter_mod (optional, int): Reject modification information
///     on any base whose basecalling quality is below this number. Defaults to 0.
/// mod_region (optional, str): Genomic region in the format "contig",
///     "contig:start-" or "contig:start-end". Reject any modification information
///     outside this region. These are half-open, 0-based intervals.
///     Can be used in combination with `region`.
///
/// # Returns
///
/// A stream of bytes that can be decoded to JSON (See the snippet from `Example output`).
///
/// # Example output
///
/// You've to decode the output of the function using something like:
///
/// ```python
/// import json
/// # Assume the function output is in 'result_bytes'
/// decoded_output = json.loads(result_bytes)
/// ```
/// A record from the decoded output might look like
///
/// ```json
/// [
/// {
///    "read_id": "cd623d4a-510d-4c6c-9d88-10eb475ac59d",
///    "sequence_length": 2104,
///    "contig": "contig_0",
///    "reference_start": 7369,
///    "reference_end": 9473,
///    "alignment_length": 2104,
///    "alignment_type": "primary_reverse",
///    "mod_count": "C-m:263;N+N:2104;(probabilities >= 0.5020, PHRED base qual >= 0)"
/// }
/// ]
/// ```
///
/// When mods are not available, you will see `NA` in the `mod_count` field.
///
/// # Errors
/// If building of option-related structs fails, if BAM input
/// cannot be obtained, if preparing records fails, or running the
/// `nanalogue_core::read_info::run` function fails
#[expect(
    clippy::doc_markdown,
    reason = "Python bindings use Python-style documentation conventions"
)]
#[expect(
    clippy::too_many_arguments,
    clippy::fn_params_excessive_bools,
    reason = "See module-level doc 'Design Philosophy: Many Parameters vs Parameter Objects'"
)]
#[pyfunction]
#[pyo3(signature = (bam_path, treat_as_url = false, min_seq_len = 0, min_align_len = 0,
                    read_ids = HashSet::<String>::new(), threads = 2, include_zero_len = false,
                    read_filter = "", sample_fraction = 1.0, mapq_filter = 0,
                    exclude_mapq_unavail = false, region = "", full_region = false,
                    tag = "", mod_strand = "", min_mod_qual = 0,
                    reject_mod_qual_non_inclusive = (0, 0), trim_read_ends_mod = 0,
                    base_qual_filter_mod = 0, mod_region = ""))]
fn read_info(
    bam_path: &str,
    treat_as_url: bool,
    min_seq_len: u64,
    min_align_len: i64,
    read_ids: HashSet<String>,
    threads: u8,
    include_zero_len: bool,
    read_filter: &str,
    sample_fraction: f32,
    mapq_filter: u8,
    exclude_mapq_unavail: bool,
    region: &str,
    full_region: bool,
    tag: &str,
    mod_strand: &str,
    min_mod_qual: u8,
    reject_mod_qual_non_inclusive: (u8, u8),
    trim_read_ends_mod: usize,
    base_qual_filter_mod: u8,
    mod_region: &str,
) -> PyResult<Vec<u8>> {
    // get input options
    let (mut bam, mut mods) = parse_input_options(
        bam_path,
        treat_as_url,
        min_seq_len,
        min_align_len,
        read_ids,
        threads,
        include_zero_len,
        read_filter,
        sample_fraction,
        mapq_filter,
        exclude_mapq_unavail,
        region,
        full_region,
        tag,
        mod_strand,
        min_mod_qual,
        reject_mod_qual_non_inclusive,
        trim_read_ends_mod,
        base_qual_filter_mod,
        mod_region,
    )?;

    // set up output buffer
    let mut buffer = Vec::new();

    // get input data and process
    let mut reader = load_bam(bam.clone())?;

    let bam_rc_records =
        BamRcRecords::new(&mut reader, &mut bam, &mut mods).map_err(|e| py_exception!(e))?;

    rust_read_info::run(
        &mut buffer,
        bam_rc_records
            .rc_records
            .filter(|r| r.as_ref().map_or(true, |v| v.pre_filt(&bam))),
        mods,
        None,
    )
    .map_err(|e| py_exception!(e))?;

    // return output
    Ok(buffer)
}

/// Window modification density across single molecules and return a Polars `DataFrame`.
/// With the gradient option, the gradient in mod density within each window is reported.
/// The output is a BED format, with windowed densities per window per read.
/// Details: runs `nanalogue_core::window_reads::run_df`.
///
/// Sets various options through builder functions before running
/// the function and capturing the output (see `Example Output`).
///
/// # Args
/// bam_path (str): Path to the BAM file. Must be associated with a BAM index.
/// win (int): Size of window in number of bases whose mod is being queried.
///     i.e. let's say a read contains cytosine mods and win is set to 10,
///     then each window is chosen so that there are 10 cytosines in it.
///     If a read has multiple mods, then multiple windows are set up such that
///     each window has the specified number of bases of that type in it.
/// step (int): Length by which the window is slid in the same units as win above.
/// win_op (optional, str): Type of windowing operation to use, allows "density" and
///     "grad_density" i.e. measure modification density or the gradient of it within
///     each window. Default is "density".
/// treat_as_url (optional, bool): If True, treat `bam_path` as a URL, default False.
/// min_seq_len (optional, int): Only retain sequences above this length, default 0.
/// min_align_len (optional, int): Only retain sequences with an alignment length above this
///     value. Defaults to unused.
/// read_ids (optional, set of str): Only retrieve these read ids, defaults to unused.
/// threads (optional, int): Number of threads used in some aspects of program execution, defaults to 2.
/// include_zero_len (optional, bool): Include sequences of zero length. WARNING: our program
///     may crash if you do this. Defaults to False. Helps to check if sequences of zero length
///     exist in our BAM file.
/// read_filter (optional, str): Comma-separated sequence of one to many of the following
///     strings: primary_forward, primary_reverse, secondary_forward, secondary_reverse,
///     supplementary_forward, supplementary_reverse, unmapped. If specified, only reads
///     with a mapping belonging to this set are retained. Defaults to no filter.
/// sample_fraction (optional, float): Set to between 0 and 1 to subsample BAM file.
///     WARNING: seeds are not set, so you may get a new set of reads every time.
///     WARNING: we sample every read with the given probability, so the total number
///         of reads fluctuates according to standard counting statistics.
/// mapq_filter (optional, int): Exclude reads with mapping quality below this number.
///     defaults to unused.
/// exclude_mapq_unavail (optional, bool): Exclude reads where mapping quality is unavailable.
///     defaults to false.
/// region (optional, str): Only include reads with at least one mapped base from this region.
///     Use the format "contig", "contig:start-", or "contig:start-end". These are 0-based,
///     half-open intervals. Defaults to read entire BAM file. Can be used in combination
///     with `mod_region`.
/// full_region (optional, bool): Only include reads if they pass through the region above
///     in full. Defaults to false.
/// tag (optional, str): If set, only this type of modification is processed. Input is a
///     string type, for example a single-letter code "m", a number as a string "76792" etc.
///     Defaults to processing all modifications.
/// mod_strand (optional, str): Set this to `bc` or `bc_comp` to retrieve information
///     about mods only from the basecalled strand or only from its complement.
///     Some sequencing technologies like `PacBio` or `ONT duplex` record mod information
///     both on a strand and its complement. It may be useful in some scenarios to
///     separate this information. Defaults to not filter.
/// min_mod_qual (optional, int): Set to a number 0-255. Reject modification
///     calls whose probability is below this value (0, 255 correspond to a
///     probability of 0 and 1 respectively). Defaults to 0.
/// reject_mod_qual_non_inclusive (optional, (int, int)): Reject modification
///     calls whose probability is such that int_low < prob < int_high.
///     Set both to a number between 0-255 and such that the first entry is <=
///     the second (if they are equal, no filtering is performed). Defaults
///     to no filtering. Also see comments under `min_mod_qual`.
/// trim_read_ends_mod (optional, int): Reject modification information
///     within so many bp of either end of the read. Defaults to 0.
/// base_qual_filter_mod (optional, int): Reject modification information
///     on any base whose basecalling quality is below this number. Defaults to 0.
/// mod_region (optional, str): Genomic region in the format "contig",
///     "contig:start-" or "contig:start-end". Reject any modification information
///     outside this region. These are half-open, 0-based intervals.
///     Can be used in combination with `region`.
///
/// # Returns
///
/// A Polars dataframe (See the snippet from `Example output`).
///
/// # Example output
///
/// If the dataframe were converted to a TSV format, it might look like the following.
/// The basecalling qualities are all 255 i.e. unknown because this is a BAM file where basecalling
/// qualities haven't been recorded.
///
/// ```text
/// #contig	ref_win_start	ref_win_end	read_id	win_val	strand	base	mod_strand	mod_type	win_start	win_end	basecall_qual
/// dummyIII	26	32	a4f36092-b4d5-47a9-813e-c22c3b477a0c	1	+	T	+	T	3	9	255
/// dummyIII	31	51	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.5	+	T	+	T	8	28	255
/// dummyIII	62	71	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.5	+	T	+	T	39	48	255
/// dummyII	22	24	fffffff1-10d2-49cb-8ca3-e8d48979001b	0.5	-	T	+	T	19	21	255
/// .	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	1	.	T	+	T	3	9	255
/// .	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.5	.	T	+	T	8	28	255
/// .	-1	-1	a4f36092-b4d5-47a9-813e-c22c3b477a0c	0.5	.	T	+	T	39	48	255
/// ```
///
/// Under the gradient option, `win_val` reports the gradient in modification density within each window.
///
/// ```text
/// #contig ref_win_start   ref_win_end     read_id win_val strand  base    mod_strand      mod_type       win_start        win_end basecall_qual
/// dummyIII    40  50  c4f36092-b4d5-47a9-813e-c22c3b477a0c    0.054545455 +   N   +   N   17  27  255
/// dummyIII    41  51  c4f36092-b4d5-47a9-813e-c22c3b477a0c    0.096969694 +   N   +   N   18  28  255
/// dummyIII    42  52  c4f36092-b4d5-47a9-813e-c22c3b477a0c    0.12727273  +   N   +   N   19  29  255
/// dummyIII    43  53  c4f36092-b4d5-47a9-813e-c22c3b477a0c    0.14545454  +   N   +   N   20  30  255
/// dummyIII    44  54  c4f36092-b4d5-47a9-813e-c22c3b477a0c    0.15151516  +   N   +   N   21  31  255
/// dummyIII    45  55  c4f36092-b4d5-47a9-813e-c22c3b477a0c    0.14545454  +   N   +   N   22  32  255
/// dummyIII    46  56  c4f36092-b4d5-47a9-813e-c22c3b477a0c    0.12727273  +   N   +   N   23  33  255
/// dummyIII    47  57  c4f36092-b4d5-47a9-813e-c22c3b477a0c    0.096969694 +   N   +   N   24  34  255
/// dummyIII    48  58  c4f36092-b4d5-47a9-813e-c22c3b477a0c    0.054545455 +   N   +   N   25  35  255
/// ```
///
/// # Errors
/// If building of option-related structs fails, if BAM input
/// cannot be obtained, if preparing records fails, or running the
/// `nanalogue_core::window_reads::run_df` function fails
#[expect(
    clippy::doc_markdown,
    reason = "Python bindings use Python-style documentation conventions"
)]
#[expect(
    clippy::too_many_arguments,
    clippy::fn_params_excessive_bools,
    reason = "See module-level doc 'Design Philosophy: Many Parameters vs Parameter Objects'"
)]
#[pyfunction]
#[pyo3(signature = (bam_path, win, step, win_op = "density", treat_as_url = false, min_seq_len = 0, min_align_len = 0,
                    read_ids = HashSet::<String>::new(), threads = 2, include_zero_len = false,
                    read_filter = "", sample_fraction = 1.0, mapq_filter = 0,
                    exclude_mapq_unavail = false, region = "", full_region = false,
                    tag = "", mod_strand = "", min_mod_qual = 0,
                    reject_mod_qual_non_inclusive = (0, 0), trim_read_ends_mod = 0,
                    base_qual_filter_mod = 0, mod_region = ""))]
fn window_reads(
    bam_path: &str,
    win: usize,
    step: usize,
    win_op: &str,
    treat_as_url: bool,
    min_seq_len: u64,
    min_align_len: i64,
    read_ids: HashSet<String>,
    threads: u8,
    include_zero_len: bool,
    read_filter: &str,
    sample_fraction: f32,
    mapq_filter: u8,
    exclude_mapq_unavail: bool,
    region: &str,
    full_region: bool,
    tag: &str,
    mod_strand: &str,
    min_mod_qual: u8,
    reject_mod_qual_non_inclusive: (u8, u8),
    trim_read_ends_mod: usize,
    base_qual_filter_mod: u8,
    mod_region: &str,
) -> PyResult<PyDataFrame> {
    // get input options
    let (mut bam, mut mods) = parse_input_options(
        bam_path,
        treat_as_url,
        min_seq_len,
        min_align_len,
        read_ids,
        threads,
        include_zero_len,
        read_filter,
        sample_fraction,
        mapq_filter,
        exclude_mapq_unavail,
        region,
        full_region,
        tag,
        mod_strand,
        min_mod_qual,
        reject_mod_qual_non_inclusive,
        trim_read_ends_mod,
        base_qual_filter_mod,
        mod_region,
    )?;

    // set up windowing options
    let win_options = InputWindowingBuilder::default()
        .win(win)
        .step(step)
        .build()
        .map_err(|e| py_value_error!(e))?;

    // get input data and process
    let mut reader = load_bam(bam.clone())?;

    let bam_rc_records =
        BamRcRecords::new(&mut reader, &mut bam, &mut mods).map_err(|e| py_exception!(e))?;

    let df = match win_op {
        "density" => rust_window_reads::run_df(
            bam_rc_records
                .rc_records
                .filter(|r| r.as_ref().map_or(true, |v| v.pre_filt(&bam))),
            win_options,
            &mods,
            |x| analysis::threshold_and_mean(x).map(Into::into),
        ),
        "grad_density" => rust_window_reads::run_df(
            bam_rc_records
                .rc_records
                .filter(|r| r.as_ref().map_or(true, |v| v.pre_filt(&bam))),
            win_options,
            &mods,
            analysis::threshold_and_gradient,
        ),
        _ => Err(Error::InvalidState(String::from(
            "win_op must be set to density or grad_density",
        ))),
    }
    .map_err(|e| py_exception!(e))?;

    // return output
    Ok(PyDataFrame(df))
}

/// Converts modification data in mod BAM files into a Polars `Dataframe`.
/// Columns are `read_id`, `seq_len`, `alignment_type`, `align_start`, `align_end`,
/// `contig`, `contig_id`, `base`, `is_strand_plus`, `mod_code`, `position`, `ref_position`,
/// and `mod_quality`.
///
/// Sets various options through builder functions before running
/// the function and capturing the output.
/// Parses records, collects them, and runs `nanalogue_core::read_utils::curr_reads_to_dataframe`.
///
/// # Args
/// bam_path (str): Path to the BAM file. Must be associated with a BAM index.
/// treat_as_url (optional, bool): If True, treat `bam_path` as a URL, default False.
/// min_seq_len (optional, int): Only retain sequences above this length, default 0.
/// min_align_len (optional, int): Only retain sequences with an alignment length above this
///     value. Defaults to unused.
/// read_ids (optional, set of str): Only retrieve these read ids, defaults to unused.
/// threads (optional, int): Number of threads used in some aspects of program execution, defaults to 2.
/// include_zero_len (optional, bool): Include sequences of zero length. WARNING: our program
///     may crash if you do this. Defaults to False. Helps to check if sequences of zero length
///     exist in our BAM file.
/// read_filter (optional, str): Comma-separated sequence of one to many of the following
///     strings: primary_forward, primary_reverse, secondary_forward, secondary_reverse,
///     supplementary_forward, supplementary_reverse, unmapped. If specified, only reads
///     with a mapping belonging to this set are retained. Defaults to no filter.
/// sample_fraction (optional, float): Set to between 0 and 1 to subsample BAM file.
///     WARNING: seeds are not set, so you may get a new set of reads every time.
///     WARNING: we sample every read with the given probability, so the total number
///         of reads fluctuates according to standard counting statistics.
/// mapq_filter (optional, int): Exclude reads with mapping quality below this number.
///     defaults to unused.
/// exclude_mapq_unavail (optional, bool): Exclude reads where mapping quality is unavailable.
///     defaults to false.
/// region (optional, str): Only include reads with at least one mapped base from this region.
///     Use the format "contig", "contig:start-", or "contig:start-end". These are 0-based,
///     half-open intervals. Defaults to read entire BAM file. Can be used in combination
///     with `mod_region`.
/// full_region (optional, bool): Only include reads if they pass through the region above
///     in full. Defaults to false.
/// tag (optional, str): If set, only this type of modification is processed. Input is a
///     string type, for example a single-letter code "m", a number as a string "76792" etc.
///     Defaults to processing all modifications.
/// mod_strand (optional, str): Set this to `bc` or `bc_comp` to retrieve information
///     about mods only from the basecalled strand or only from its complement.
///     Some sequencing technologies like `PacBio` or `ONT duplex` record mod information
///     both on a strand and its complement. It may be useful in some scenarios to
///     separate this information. Defaults to not filter.
/// min_mod_qual (optional, int): Set to a number 0-255. Reject modification
///     calls whose probability is below this value (0, 255 correspond to a
///     probability of 0 and 1 respectively). Defaults to 0.
/// reject_mod_qual_non_inclusive (optional, (int, int)): Reject modification
///     calls whose probability is such that int_low < prob < int_high.
///     Set both to a number between 0-255 and such that the first entry is <=
///     the second (if they are equal, no filtering is performed). Defaults
///     to no filtering. Also see comments under `min_mod_qual`.
/// trim_read_ends_mod (optional, int): Reject modification information
///     within so many bp of either end of the read. Defaults to 0.
/// base_qual_filter_mod (optional, int): Reject modification information
///     on any base whose basecalling quality is below this number. Defaults to 0.
/// mod_region (optional, str): Genomic region in the format "contig",
///     "contig:start-" or "contig:start-end". Reject any modification information
///     outside this region. These are half-open, 0-based intervals.
///     Can be used in combination with `region`.
///
/// # Returns
///
/// A Polars dataframe.
///
/// # Example output
///
/// If the polars dataframe were converted to a plain-text table, it might look like the following:
///
/// ```text
/// read_id	seq_len	alignment_type	align_start	align_end	contig	contig_id	base	is_strand_plus	mod_code	position	ref_position	mod_quality
/// 5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	8	primary_forward	9	17	dummyI	0	T	true	T	0	9	4
/// 5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	8	primary_forward	9	17	dummyI	0	T	true	T	3	12	7
/// 5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	8	primary_forward	9	17	dummyI	0	T	true	T	4	13	9
/// 5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	8	primary_forward	9	17	dummyI	0	T	true	T	7	16	6
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	primary_forward	23	71	dummyIII	2	T	true	T	3	26	221
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	primary_forward	23	71	dummyIII	2	T	true	T	8	31	242
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	primary_forward	23	71	dummyIII	2	T	true	T	27	50	3
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	primary_forward	23	71	dummyIII	2	T	true	T	39	62	47
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	primary_forward	23	71	dummyIII	2	T	true	T	47	70	239
/// fffffff1-10d2-49cb-8ca3-e8d48979001b	33	primary_reverse	3	36	dummyII	1	T	true	T	12	15	3
/// fffffff1-10d2-49cb-8ca3-e8d48979001b	33	primary_reverse	3	36	dummyII	1	T	true	T	13	16	3
/// fffffff1-10d2-49cb-8ca3-e8d48979001b	33	primary_reverse	3	36	dummyII	1	T	true	T	16	19	4
/// fffffff1-10d2-49cb-8ca3-e8d48979001b	33	primary_reverse	3	36	dummyII	1	T	true	T	19	22	3
/// fffffff1-10d2-49cb-8ca3-e8d48979001b	33	primary_reverse	3	36	dummyII	1	T	true	T	20	23	182
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					G	false	7200	28	-1	0
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					G	false	7200	29	-1	0
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					G	false	7200	30	-1	0
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					G	false	7200	32	-1	0
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					G	false	7200	43	-1	77
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					G	false	7200	44	-1	0
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					T	true	T	3	-1	221
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					T	true	T	8	-1	242
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					T	true	T	27	-1	0
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					T	true	T	39	-1	47
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	48	unmapped					T	true	T	47	-1	239
/// ```
///
/// # Errors
/// If building of option-related structs fails, if BAM input
/// cannot be obtained, if preparing records fails, or running the
/// `nanalogue_core` commands fail
#[expect(
    clippy::doc_markdown,
    reason = "Python bindings use Python-style documentation conventions"
)]
#[expect(
    clippy::too_many_arguments,
    clippy::fn_params_excessive_bools,
    reason = "See module-level doc 'Design Philosophy: Many Parameters vs Parameter Objects'"
)]
#[pyfunction]
#[pyo3(signature = (bam_path, treat_as_url = false, min_seq_len = 0, min_align_len = 0,
                    read_ids = HashSet::<String>::new(), threads = 2, include_zero_len = false,
                    read_filter = "", sample_fraction = 1.0, mapq_filter = 0,
                    exclude_mapq_unavail = false, region = "", full_region = false,
                    tag = "", mod_strand = "", min_mod_qual = 0,
                    reject_mod_qual_non_inclusive = (0, 0), trim_read_ends_mod = 0,
                    base_qual_filter_mod = 0, mod_region = ""))]
fn polars_bam_mods(
    bam_path: &str,
    treat_as_url: bool,
    min_seq_len: u64,
    min_align_len: i64,
    read_ids: HashSet<String>,
    threads: u8,
    include_zero_len: bool,
    read_filter: &str,
    sample_fraction: f32,
    mapq_filter: u8,
    exclude_mapq_unavail: bool,
    region: &str,
    full_region: bool,
    tag: &str,
    mod_strand: &str,
    min_mod_qual: u8,
    reject_mod_qual_non_inclusive: (u8, u8),
    trim_read_ends_mod: usize,
    base_qual_filter_mod: u8,
    mod_region: &str,
) -> PyResult<PyDataFrame> {
    // get input options
    let (mut bam, mut mods) = parse_input_options(
        bam_path,
        treat_as_url,
        min_seq_len,
        min_align_len,
        read_ids,
        threads,
        include_zero_len,
        read_filter,
        sample_fraction,
        mapq_filter,
        exclude_mapq_unavail,
        region,
        full_region,
        tag,
        mod_strand,
        min_mod_qual,
        reject_mod_qual_non_inclusive,
        trim_read_ends_mod,
        base_qual_filter_mod,
        mod_region,
    )?;

    // prepare output
    let mut df_collection = Vec::new();

    // get input data and process
    let mut reader = load_bam(bam.clone())?;

    let bam_rc_records =
        BamRcRecords::new(&mut reader, &mut bam, &mut mods).map_err(|e| py_exception!(e))?;

    for k in bam_rc_records
        .rc_records
        .filter(|r| r.as_ref().map_or(true, |v| v.pre_filt(&bam)))
    {
        let record = k.map_err(|e| py_exception!(e))?;
        let curr_read = CurrRead::default()
            .try_from_only_alignment(&record)
            .map_err(|e| py_exception!(e))?
            .set_mod_data_restricted_options(&record, &mods)
            .map_err(|e| py_exception!(e))?;
        df_collection.push(curr_read);
    }

    // Convert to DataFrame
    let df =
        curr_reads_to_dataframe(df_collection.as_slice()).map_err(|e: Error| py_exception!(e))?;

    // return output
    Ok(PyDataFrame(df))
}

/// Returns a Polars `DataFrame` with read sequence and quality information for a genomic region.
/// Represents insertions as lowercase, deletions as periods, and modifications as 'Z'.
/// Basecalling qualities represented as a string of numbers separated by periods, with value set
/// to '255' at a deleted base.
///
/// Calls `nanalogue_core::reads_table::run_df` to extract read sequences and qualities
/// from a specified genomic region. The output DataFrame contains columns for read ID,
/// sequence (with modifications shown as 'Z'), and base qualities.
///
/// # Args
/// bam_path (str): Path to the BAM file. Must be associated with a BAM index.
/// region (str): Genomic region from which to extract sequences. Required.
///     Use the format "contig", "contig:start-", or "contig:start-end". These are 0-based,
///     half-open intervals.
/// treat_as_url (optional, bool): If True, treat `bam_path` as a URL, default False.
/// min_seq_len (optional, int): Only retain sequences above this length, default 0.
/// min_align_len (optional, int): Only retain sequences with an alignment length above this
///     value. Defaults to unused.
/// read_ids (optional, set of str): Only retrieve these read ids, defaults to unused.
/// threads (optional, int): Number of threads used in some aspects of program execution, defaults to 2.
/// include_zero_len (optional, bool): Include sequences of zero length. WARNING: our program
///     may crash if you do this. Defaults to False.
/// read_filter (optional, str): Comma-separated sequence of one to many of the following
///     strings: primary_forward, primary_reverse, secondary_forward, secondary_reverse,
///     supplementary_forward, supplementary_reverse, unmapped. If specified, only reads
///     with a mapping belonging to this set are retained. Defaults to no filter.
/// sample_fraction (optional, float): Set to between 0 and 1 to subsample BAM file.
///     WARNING: seeds are not set, so you may get a new set of reads every time.
/// mapq_filter (optional, int): Exclude reads with mapping quality below this number.
///     defaults to unused.
/// exclude_mapq_unavail (optional, bool): Exclude reads where mapping quality is unavailable.
///     defaults to false.
/// tag (optional, str): If set, only this type of modification is processed. Input is a
///     string type, for example a single-letter code "m", a number as a string "76792" etc.
///     Defaults to processing all modifications.
/// mod_strand (optional, str): Set this to `bc` or `bc_comp` to retrieve information
///     about mods only from the basecalled strand or only from its complement.
///     Defaults to not filter.
/// min_mod_qual (optional, int): Set to a number 0-255. Reject modification
///     calls whose probability is below this value. Defaults to 0.
/// reject_mod_qual_non_inclusive (optional, (int, int)): Reject modification
///     calls whose probability is such that int_low < prob < int_high.
///     Defaults to no filtering.
/// trim_read_ends_mod (optional, int): Reject modification information
///     within so many bp of either end of the read. Defaults to 0.
/// base_qual_filter_mod (optional, int): Reject modification information
///     on any base whose basecalling quality is below this number. Defaults to 0.
///
/// # Returns
///
/// A Polars dataframe with columns: `read_id`, `sequence`, `qualities`.
///
/// Sequence column conventions:
/// - Insertions are shown in lowercase
/// - Deletions are shown as periods (`.`)
/// - Modified bases on the reference are shown as `Z`
/// - Modified bases in an insertion are shown as `z`
///
/// Qualities column: a period-separated string of integer quality scores,
/// with one score per base in the sequence.
///
/// # Example output
///
/// If the polars dataframe were converted to a plain-text table, it might look like:
///
/// ```text
/// read_id	sequence	qualities
/// 5d10eb9a-aae1-4db8-8ec6-7ebb34d32575	ACGTZac..TZgt	30.35.40.38.42.20.22.255.255.41.25.28.30
/// a4f36092-b4d5-47a9-813e-c22c3b477a0c	TGCAZz.ATGCA	28.33.39.41.40.18.255.36.42.38.35.30
/// ```
///
/// # Errors
/// If the region parameter is empty, if building of option-related structs fails,
/// if BAM input cannot be obtained, if preparing records fails, or if running the
/// `nanalogue_core::reads_table::run_df` function fails
#[expect(
    clippy::doc_markdown,
    reason = "Python bindings use Python-style documentation conventions"
)]
#[expect(
    clippy::too_many_arguments,
    reason = "See module-level doc 'Design Philosophy: Many Parameters vs Parameter Objects'"
)]
#[pyfunction]
#[pyo3(signature = (bam_path, region, treat_as_url = false, min_seq_len = 0, min_align_len = 0,
                    read_ids = HashSet::<String>::new(), threads = 2, include_zero_len = false,
                    read_filter = "", sample_fraction = 1.0, mapq_filter = 0,
                    exclude_mapq_unavail = false, tag = "", mod_strand = "", min_mod_qual = 0,
                    reject_mod_qual_non_inclusive = (0, 0), trim_read_ends_mod = 0,
                    base_qual_filter_mod = 0))]
fn seq_table(
    bam_path: &str,
    region: &str,
    treat_as_url: bool,
    min_seq_len: u64,
    min_align_len: i64,
    read_ids: HashSet<String>,
    threads: u8,
    include_zero_len: bool,
    read_filter: &str,
    sample_fraction: f32,
    mapq_filter: u8,
    exclude_mapq_unavail: bool,
    tag: &str,
    mod_strand: &str,
    min_mod_qual: u8,
    reject_mod_qual_non_inclusive: (u8, u8),
    trim_read_ends_mod: usize,
    base_qual_filter_mod: u8,
) -> PyResult<PyDataFrame> {
    // region is required for SeqDisplayOptions::Region
    if region.is_empty() {
        return Err(py_value_error!(
            "region parameter is required for seq_table (cannot be empty)"
        ));
    }

    // get input options with full_region = true and mod_region = region
    let (mut bam, mut mods) = parse_input_options(
        bam_path,
        treat_as_url,
        min_seq_len,
        min_align_len,
        read_ids,
        threads,
        include_zero_len,
        read_filter,
        sample_fraction,
        mapq_filter,
        exclude_mapq_unavail,
        region,
        true, // full_region hardcoded to true
        tag,
        mod_strand,
        min_mod_qual,
        reject_mod_qual_non_inclusive,
        trim_read_ends_mod,
        base_qual_filter_mod,
        region, // mod_region same as region
    )?;

    // get input data
    let mut reader = load_bam(bam.clone())?;

    let bam_rc_records =
        BamRcRecords::new(&mut reader, &mut bam, &mut mods).map_err(|e| py_exception!(e))?;

    // Parse region to GenomicRegion then convert to Bed3 for SeqDisplayOptions
    let genomic_region = GenomicRegion::from_str(region).map_err(|e: Error| py_value_error!(e))?;
    let region_bed3 = genomic_region
        .try_to_bed3(&bam_rc_records.header)
        .map_err(|e: Error| py_value_error!(e))?;

    // Create SeqDisplayOptions::Region with specified settings
    let seq_display = SeqDisplayOptions::Region {
        show_base_qual: true,
        show_ins_lowercase: true,
        region: region_bed3,
        show_mod_z: true,
    };

    // Call reads_table::run_df with seq_summ_path = ""
    let full_df = rust_reads_table::run_df(
        bam_rc_records
            .rc_records
            .filter(|r| r.as_ref().map_or(true, |v| v.pre_filt(&bam))),
        Some(mods),
        seq_display,
        "", // seq_summ_path empty
    )
    .map_err(|e| py_exception!(e))?;

    // Select only the required columns: read_id, sequence, qualities
    let df = full_df
        .select(["read_id", "sequence", "qualities"])
        .map_err(|e| py_exception!(e))?;

    // return output
    Ok(PyDataFrame(df))
}

/// Simulates a BAM file with or without modifications based on a JSON configuration.
/// Creates both a BAM file and a corresponding FASTA reference file.
///
/// This function takes a JSON configuration string that specifies how to generate
/// synthetic sequencing data, including contig specifications, read parameters,
/// and optional modification patterns. The output includes a sorted, indexed BAM
/// file and a FASTA reference file.
///
/// # Args
/// json_config (str): JSON string containing the simulation configuration.
///     Must conform to the `SimulationConfig` schema. The JSON should specify
///     (some inputs may be optional):
///     - `contigs`: Configuration for generating reference sequences
///         - `number`: Number of contigs to generate
///         - `len_range`: Range of contig lengths [min, max]
///         - `repeated_seq`: Sequence pattern to repeat for contig generation.
///             If not specified, random DNA sequences are generated.
///     - `reads`: Array of read group specifications, each containing:
///         - `number`: Number of reads to generate
///         - `mapq_range`: Range of mapping quality values [min, max]
///         - `base_qual_range`: Range of base quality values [min, max]
///         - `len_range`: Range of read lengths as fraction of contig [min, max]
///         - `barcode`: Optional barcode sequence
///         - `mods`: Optional array of modification specifications
/// bam_path (str): Output path for the BAM file. If the file exists, it will be overwritten.
///     A corresponding `.bai` index file will be created automatically.
/// fasta_path (str): Output path for the FASTA reference file. If the file exists,
///     it will be overwritten.
///
/// # Returns
///
/// Returns None on success. The function creates two files on disk:
/// - A sorted, indexed BAM file at `bam_path` (with accompanying `.bai` index)
/// - A FASTA reference file at `fasta_path`
///
/// # Example
///
/// ```python
/// import pynanalogue
///
/// json_config = '''
/// {
/// "contigs": {
///     "number": 2,
///     "len_range": [100, 200],
///     "repeated_seq": "ACGTACGT"
/// },
/// "reads": [
///     {
///         "number": 10,
///         "mapq_range": [10, 30],
///         "base_qual_range": [20, 40],
///         "len_range": [0.1, 0.9],
///         "barcode": "ACGTAA",
///         "mods": [{
///             "base": "C",
///             "is_strand_plus": true,
///             "mod_code": "m",
///             "win": [5, 3],
///             "mod_range": [[0.3, 0.7], [0.1, 0.5]]
///         }]
///     }
/// ]
/// }
/// '''
///
/// pynanalogue.simulate_mod_bam(
/// json_config=json_config,
/// bam_path="output.bam",
/// fasta_path="output.fasta"
/// )
/// ```
///
/// # Errors
/// Returns a Python exception if:
/// - The JSON configuration is invalid or cannot be parsed
/// - The JSON structure doesn't match the expected `SimulationConfig` schema
/// - File I/O operations fail (e.g., permission issues, disk full)
/// - BAM or FASTA generation fails due to invalid configuration parameters
#[expect(
    clippy::doc_markdown,
    reason = "Python bindings use Python-style documentation conventions"
)]
#[pyfunction]
#[pyo3(signature = (json_config, bam_path, fasta_path))]
fn simulate_mod_bam(json_config: &str, bam_path: &str, fasta_path: &str) -> PyResult<()> {
    // Parse JSON string into SimulationConfig
    let config: SimulationConfig =
        serde_json::from_str(json_config).map_err(|e| py_value_error!(e))?;

    // Run the simulation
    rust_simulate_mod_bam::run(config, bam_path, fasta_path).map_err(|e| py_exception!(e))?;

    Ok(())
}

/// Peeks at a BAM file to extract contig information and detected modifications.
///
/// This function reads the BAM header and examines up to 100 records to determine
/// the contigs present in the file and any DNA/RNA modifications detected.
///
/// # Args
/// bam_path (str): Path to the BAM file.
/// treat_as_url (optional, bool): If True, treat `bam_path` as a URL, default False.
///
/// # Returns
///
/// A dictionary with two keys:
/// - `contigs`: A dictionary mapping contig names to their lengths
/// - `modifications`: A list of modifications, where each modification is a list
///   of three strings: \[base, strand, `mod_code`\]. For example: \[\["T", "+", "T"\], \["G", "-", "7200"\]\]
///
/// # Example
///
/// ```python
/// import pynanalogue
///
/// result = pynanalogue.peek("example.bam")
/// print(result["contigs"])  # {"chr1": 248956422, "chr2": 242193529, ...}
/// print(result["modifications"])  # [["C", "+", "m"], ["A", "+", "28871"]]
/// ```
///
/// # Errors
/// If the BAM file cannot be read or parsed
#[expect(
    clippy::doc_markdown,
    reason = "Python bindings use Python-style documentation conventions"
)]
#[pyfunction]
#[pyo3(signature = (bam_path, treat_as_url = false))]
fn peek(py: Python<'_>, bam_path: &str, treat_as_url: bool) -> PyResult<Py<PyDict>> {
    // Build minimal InputBam with just path/url
    let mut input_bam = InputBamBuilder::default()
        .bam_path({
            if treat_as_url {
                PathOrURLOrStdin::URL(
                    Url::parse(bam_path).map_err(|e: ParseError| py_value_error!(e))?,
                )
            } else {
                PathOrURLOrStdin::Path(bam_path.into())
            }
        })
        .build()
        .map_err(|e| py_value_error!(e))?;

    // Load BAM and create BamRcRecords to get header and records
    let mut reader = load_bam(input_bam.clone())?;
    let bam_rc_records = BamRcRecords::new(
        &mut reader,
        &mut input_bam,
        &mut InputMods::<OptionalTag>::default(),
    )
    .map_err(|e| py_exception!(e))?;

    // Run peek and capture output
    let mut buffer = Vec::new();
    rust_peek::run(
        &mut buffer,
        &bam_rc_records.header,
        bam_rc_records.rc_records.take(100),
    )
    .map_err(|e| py_exception!(e))?;

    // Parse the output string directly into Python dict/list
    let output_str = String::from_utf8(buffer).map_err(|e| py_value_error!(e))?;
    let result = PyDict::new(py);
    let contigs_dict = PyDict::new(py);
    let mods_list = PyList::empty(py);
    let mut in_contigs_section = true;

    for line in output_str.lines() {
        let trimmed = line.trim();
        match trimmed {
            "" | "None" => {}
            "contigs_and_lengths:" => in_contigs_section = true,
            "modifications:" => in_contigs_section = false,
            _ if in_contigs_section => {
                // Parse "contig_name\tlength"
                let parts: Vec<&str> = trimmed.split('\t').collect();
                let contig_name = parts
                    .first()
                    .ok_or_else(|| py_value_error!("Missing contig name"))?;
                let length: u64 = parts
                    .get(1)
                    .ok_or_else(|| py_value_error!("Missing contig length"))?
                    .parse()
                    .map_err(|e| py_value_error!(format!("Failed to parse contig length: {e}")))?;
                contigs_dict.set_item(*contig_name, length)?;
            }
            _ => {
                // Parse modification string like "G-7200" or "T+T"
                // Format: base + strand + mod_code (strand is always '+' or '-' at position 1)
                let base = trimmed
                    .chars()
                    .next()
                    .ok_or_else(|| py_value_error!("Empty modification string"))?
                    .to_string();
                let strand = trimmed
                    .chars()
                    .nth(1)
                    .ok_or_else(|| py_value_error!("Modification string missing strand"))?
                    .to_string();
                let mod_code: String = trimmed.chars().skip(2).collect();
                mods_list.append(PyList::new(py, [base, strand, mod_code])?)?;
            }
        }
    }

    result.set_item("contigs", contigs_dict)?;
    result.set_item("modifications", mods_list)?;

    Ok(result.into())
}

/// Our python module; calls our rust functions
///
/// # Errors
/// `PyO3` errors
#[pymodule]
fn pynanalogue(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_info, m)?)?;
    m.add_function(wrap_pyfunction!(window_reads, m)?)?;
    m.add_function(wrap_pyfunction!(polars_bam_mods, m)?)?;
    m.add_function(wrap_pyfunction!(seq_table, m)?)?;
    m.add_function(wrap_pyfunction!(simulate_mod_bam, m)?)?;
    m.add_function(wrap_pyfunction!(peek, m)?)?;

    Ok(())
}
