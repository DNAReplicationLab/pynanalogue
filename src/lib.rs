//! # Pynanalogue (Python functions for Nanalogue (= Nucleic Acid Analogue))
//!
//! In Nanalogue, we process and analyse data associated with DNA/RNA molecules,
//! their alignments to reference genomes, modification information on them,
//! and other miscellaneous information from BAM files. We expose some of
//! these functions in a python package for usage by the python bioinformatics
//! community.
use nanalogue_core::{
    BamPreFilt as _, BamRcRecords, CurrRead, Error, F32Bw0and1, InputBam, InputBamBuilder,
    InputMods, InputModsBuilder, InputWindowingBuilder, OptionalTag, OrdPair, PathOrURLOrStdin,
    ThresholdState, analysis, curr_reads_to_dataframe, nanalogue_indexed_bam_reader,
    nanalogue_indexed_bam_reader_from_url, read_info as rust_read_info,
    window_reads as rust_window_reads,
};
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use rust_htslib::bam::FetchDefinition;
use std::collections::HashSet;
use std::num::NonZeroU32;
use url::{ParseError, Url};

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
#[expect(
    clippy::let_underscore_untyped,
    reason = "occasionally we will leave _ untyped"
)]
#[expect(
    clippy::too_many_arguments,
    clippy::fn_params_excessive_bools,
    reason = "we have no choice, python allows many args (most are optional). So we've to resort to this!"
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
    mod_strand: &str,
    min_mod_qual: u8,
    reject_mod_qual_non_inclusive: (u8, u8),
    trim_read_ends_mod: usize,
    base_qual_filter_mod: u8,
    mod_region: &str,
) -> PyResult<(InputBam, InputMods<OptionalTag>)> {
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
        .read_filter(read_filter.into())
        .sample_fraction(F32Bw0and1::new(sample_fraction).map_err(|e| py_value_error!(e))?)
        .mapq_filter(mapq_filter)
        .exclude_mapq_unavail(exclude_mapq_unavail)
        .region(region.into())
        .full_region(full_region);
    let _: Option<&mut _> = (min_align_len > 0).then(|| bam_builder.min_align_len(min_align_len));
    let _: Option<&mut _> = (!read_ids.is_empty()).then(|| bam_builder.read_id_set(read_ids));

    let bam = bam_builder.build().map_err(|e| py_value_error!(e))?;

    #[expect(
        clippy::arithmetic_side_effects,
        reason = "we check low < high - 1 so no chance of over or underflow"
    )]
    let mods = InputModsBuilder::<OptionalTag>::default()
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
                _ => ThresholdState::Both((
                    min_mod_qual,
                    OrdPair::<u8>::try_from((low + 1, high - 1)).map_err(|e| py_value_error!(e))?,
                )),
            }
        })
        .trim_read_ends_mod(trim_read_ends_mod)
        .base_qual_filter_mod(base_qual_filter_mod)
        .mod_region(mod_region.into())
        .build()
        .map_err(|e| py_value_error!(e))?;

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
        _ => unreachable!("we don't allow the PathOrURLOrStdin::Stdin variant here"),
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
///     bam_path (str): Path to the BAM file. Must be associated with a BAM index.
///     treat_as_url (optional, bool): If True, treat `bam_path` as a URL, default False.
///     min_seq_len (optional, int): Only retain sequences above this length, default 0.
///     min_align_len (optional, int): Only retain sequences with an alignment length above this
///         value. Defaults to unused.
///     read_ids (optional, set of str): Only retrieve these read ids, defaults to unused.
///     threads (optional, int): Number of threads used in some aspects of program execution, defaults to 2.
///     include_zero_len (optional, bool): Include sequences of zero length. WARNING: our program
///         may crash if you do this. Defaults to False. Helps to check if sequences of zero length
///         exist in our BAM file.
///     read_filter (optional, str): Comma-separated sequence of one to many of the following
///         strings: primary_forward, primary_reverse, secondary_forward, secondary_reverse,
///         supplementary_forward, supplementary_reverse, unmapped. If specified, only reads
///         with a mapping belonging to this set are retained. Defaults to no filter.
///     sample_fraction (optional, float): Set to between 0 and 1 to subsample BAM file.
///         WARNING: seeds are not set, so you may get a new set of reads every time.
///         WARNING: we sample every read with the given probability, so the total number
///             of reads fluctuates according to standard counting statistics.
///     mapq_filter (optional, int): Exclude reads with mapping quality below this number.
///         defaults to unused.
///     exclude_mapq_unavail (optional, bool): Exclude reads where mapping quality is unavailable.
///         defaults to false.
///     region (optional, str): Only include reads with at least one mapped base from this region.
///         Use the format "contig", "contig:start-", or "contig:start-end". These are 0-based,
///         half-open intervals. Defaults to read entire BAM file. Can be used in combination
///         with `mod_region`.
///     full_region (optional, bool): Only include reads if they pass through the region above
///         in full. Defaults to false.
///     mod_strand (optional, str): Set this to `bc` or `bc_comp` to retrieve information
///         about mods only from the basecalled strand or only from its complement.
///         Some sequencing technologies like `PacBio` or `ONT duplex` record mod information
///         both on a strand and its complement. It may be useful in some scenarios to
///         separate this information. Defaults to not filter.
///     min_mod_qual (optional, int): Set to a number 0-255. Reject modification
///         calls whose probability is below this value (0, 255 correspond to a
///         probability of 0 and 1 respectively). Defaults to 0.
///     reject_mod_qual_non_inclusive (optional, (int, int)): Reject modification
///         calls whose probability is such that int_low < prob < int_high.
///         Set both to a number between 0-255 and such that the first entry is <=
///         the second (if they are equal, no filtering is performed). Defaults
///         to no filtering. Also see comments under `min_mod_qual`.
///     trim_read_ends_mod (optional, int): Reject modification information
///         within so many bp of either end of the read. Defaults to 0.
///     base_qual_filter_mod (optional, int): Reject modification information
///         on any base whose basecalling quality is below this number. Defaults to 0.
///     mod_region (optional, str): Genomic region in the format "contig",
///         "contig:start-" or "contig:start-end". Reject any modification information
///         outside this region. These are half-open, 0-based intervals.
///         Can be used in combination with `region`.
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
/// // assume the function output is in out
/// decoded_output = json.loads(bytearray(out))
/// ```
/// A record from the decoded output might look like
///
/// ```json
/// [
/// {
///        "read_id": "cd623d4a-510d-4c6c-9d88-10eb475ac59d",
///        "sequence_length": 2104,
///        "contig": "contig_0",
///        "reference_start": 7369,
///        "reference_end": 9473,
///        "alignment_length": 2104,
///        "alignment_type": "primary_reverse",
///        "mod_count": "C-m:263;N+N:2104;(probabilities >= 0.5020, PHRED base qual >= 0)"
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
    clippy::too_many_arguments,
    clippy::fn_params_excessive_bools,
    reason = "python functions have many more args than rust, o.k. as some are optional"
)]
#[pyfunction]
#[pyo3(signature = (bam_path, treat_as_url = false, min_seq_len = 0, min_align_len = 0,
                    read_ids = HashSet::<String>::new(), threads = 2, include_zero_len = false,
                    read_filter = "", sample_fraction = 1.0, mapq_filter = 0,
                    exclude_mapq_unavail = false, region = "", full_region = false,
                    mod_strand = "", min_mod_qual = 0,
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
        &mods,
        None,
    )
    .map_err(|e| py_exception!(e))?;

    // return output
    Ok(buffer)
}

/// Runs `nanalogue_core::window_reads::run_df` and gets output as a Polars `DataFrame`.
/// The output is a BED format, with windowed densities per read.
///
/// Sets various options through builder functions before running
/// the function and capturing the output (see `Example Output`).
///
/// # Args
///     bam_path (str): Path to the BAM file. Must be associated with a BAM index.
///     win (int): Size of window in number of bases whose mod is being queried.
///         i.e. let's say a read contains cytosine mods and win is set to 10,
///         then each window is chosen so that there are 10 cytosines in it.
///         If a read has multiple mods, then multiple windows are set up such that
///         each window has the specified number of bases of that type in it.
///     step (int): Length by which the window is slid in the same units as win above.
///     treat_as_url (optional, bool): If True, treat `bam_path` as a URL, default False.
///     min_seq_len (optional, int): Only retain sequences above this length, default 0.
///     min_align_len (optional, int): Only retain sequences with an alignment length above this
///         value. Defaults to unused.
///     read_ids (optional, set of str): Only retrieve these read ids, defaults to unused.
///     threads (optional, int): Number of threads used in some aspects of program execution, defaults to 2.
///     include_zero_len (optional, bool): Include sequences of zero length. WARNING: our program
///         may crash if you do this. Defaults to False. Helps to check if sequences of zero length
///         exist in our BAM file.
///     read_filter (optional, str): Comma-separated sequence of one to many of the following
///         strings: primary_forward, primary_reverse, secondary_forward, secondary_reverse,
///         supplementary_forward, supplementary_reverse, unmapped. If specified, only reads
///         with a mapping belonging to this set are retained. Defaults to no filter.
///     sample_fraction (optional, float): Set to between 0 and 1 to subsample BAM file.
///         WARNING: seeds are not set, so you may get a new set of reads every time.
///         WARNING: we sample every read with the given probability, so the total number
///             of reads fluctuates according to standard counting statistics.
///     mapq_filter (optional, int): Exclude reads with mapping quality below this number.
///         defaults to unused.
///     exclude_mapq_unavail (optional, bool): Exclude reads where mapping quality is unavailable.
///         defaults to false.
///     region (optional, str): Only include reads with at least one mapped base from this region.
///         Use the format "contig", "contig:start-", or "contig:start-end". These are 0-based,
///         half-open intervals. Defaults to read entire BAM file. Can be used in combination
///         with `mod_region`.
///     full_region (optional, bool): Only include reads if they pass through the region above
///         in full. Defaults to false.
///     mod_strand (optional, str): Set this to `bc` or `bc_comp` to retrieve information
///         about mods only from the basecalled strand or only from its complement.
///         Some sequencing technologies like `PacBio` or `ONT duplex` record mod information
///         both on a strand and its complement. It may be useful in some scenarios to
///         separate this information. Defaults to not filter.
///     min_mod_qual (optional, int): Set to a number 0-255. Reject modification
///         calls whose probability is below this value (0, 255 correspond to a
///         probability of 0 and 1 respectively). Defaults to 0.
///     reject_mod_qual_non_inclusive (optional, (int, int)): Reject modification
///         calls whose probability is such that int_low < prob < int_high.
///         Set both to a number between 0-255 and such that the first entry is <=
///         the second (if they are equal, no filtering is performed). Defaults
///         to no filtering. Also see comments under `min_mod_qual`.
///     trim_read_ends_mod (optional, int): Reject modification information
///         within so many bp of either end of the read. Defaults to 0.
///     base_qual_filter_mod (optional, int): Reject modification information
///         on any base whose basecalling quality is below this number. Defaults to 0.
///     mod_region (optional, str): Genomic region in the format "contig",
///         "contig:start-" or "contig:start-end". Reject any modification information
///         outside this region. These are half-open, 0-based intervals.
///         Can be used in combination with `region`.
///
/// # Returns
///
/// A Polars dataframe (See the snippet from `Example output`).
///
/// # Example output
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
/// # Errors
/// If building of option-related structs fails, if BAM input
/// cannot be obtained, if preparing records fails, or running the
/// `nanalogue_core::window_reads::run_df` function fails
#[expect(
    clippy::too_many_arguments,
    clippy::fn_params_excessive_bools,
    reason = "python functions have many more args than rust, o.k. as some are optional"
)]
#[pyfunction]
#[pyo3(signature = (bam_path, win, step, treat_as_url = false, min_seq_len = 0, min_align_len = 0,
                    read_ids = HashSet::<String>::new(), threads = 2, include_zero_len = false,
                    read_filter = "", sample_fraction = 1.0, mapq_filter = 0,
                    exclude_mapq_unavail = false, region = "", full_region = false,
                    mod_strand = "", min_mod_qual = 0,
                    reject_mod_qual_non_inclusive = (0, 0), trim_read_ends_mod = 0,
                    base_qual_filter_mod = 0, mod_region = ""))]
fn window_reads(
    bam_path: &str,
    win: usize,
    step: usize,
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

    let df = rust_window_reads::run_df(
        bam_rc_records
            .rc_records
            .filter(|r| r.as_ref().map_or(true, |v| v.pre_filt(&bam))),
        win_options,
        &mods,
        |x| analysis::threshold_and_mean(x).map(Into::into),
    )
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
///     bam_path (str): Path to the BAM file. Must be associated with a BAM index.
///     treat_as_url (optional, bool): If True, treat `bam_path` as a URL, default False.
///     min_seq_len (optional, int): Only retain sequences above this length, default 0.
///     min_align_len (optional, int): Only retain sequences with an alignment length above this
///         value. Defaults to unused.
///     read_ids (optional, set of str): Only retrieve these read ids, defaults to unused.
///     threads (optional, int): Number of threads used in some aspects of program execution, defaults to 2.
///     include_zero_len (optional, bool): Include sequences of zero length. WARNING: our program
///         may crash if you do this. Defaults to False. Helps to check if sequences of zero length
///         exist in our BAM file.
///     read_filter (optional, str): Comma-separated sequence of one to many of the following
///         strings: primary_forward, primary_reverse, secondary_forward, secondary_reverse,
///         supplementary_forward, supplementary_reverse, unmapped. If specified, only reads
///         with a mapping belonging to this set are retained. Defaults to no filter.
///     sample_fraction (optional, float): Set to between 0 and 1 to subsample BAM file.
///         WARNING: seeds are not set, so you may get a new set of reads every time.
///         WARNING: we sample every read with the given probability, so the total number
///             of reads fluctuates according to standard counting statistics.
///     mapq_filter (optional, int): Exclude reads with mapping quality below this number.
///         defaults to unused.
///     exclude_mapq_unavail (optional, bool): Exclude reads where mapping quality is unavailable.
///         defaults to false.
///     region (optional, str): Only include reads with at least one mapped base from this region.
///         Use the format "contig", "contig:start-", or "contig:start-end". These are 0-based,
///         half-open intervals. Defaults to read entire BAM file. Can be used in combination
///         with `mod_region`.
///     full_region (optional, bool): Only include reads if they pass through the region above
///         in full. Defaults to false.
///     mod_strand (optional, str): Set this to `bc` or `bc_comp` to retrieve information
///         about mods only from the basecalled strand or only from its complement.
///         Some sequencing technologies like `PacBio` or `ONT duplex` record mod information
///         both on a strand and its complement. It may be useful in some scenarios to
///         separate this information. Defaults to not filter.
///     min_mod_qual (optional, int): Set to a number 0-255. Reject modification
///         calls whose probability is below this value (0, 255 correspond to a
///         probability of 0 and 1 respectively). Defaults to 0.
///     reject_mod_qual_non_inclusive (optional, (int, int)): Reject modification
///         calls whose probability is such that int_low < prob < int_high.
///         Set both to a number between 0-255 and such that the first entry is <=
///         the second (if they are equal, no filtering is performed). Defaults
///         to no filtering. Also see comments under `min_mod_qual`.
///     trim_read_ends_mod (optional, int): Reject modification information
///         within so many bp of either end of the read. Defaults to 0.
///     base_qual_filter_mod (optional, int): Reject modification information
///         on any base whose basecalling quality is below this number. Defaults to 0.
///     mod_region (optional, str): Genomic region in the format "contig",
///         "contig:start-" or "contig:start-end". Reject any modification information
///         outside this region. These are half-open, 0-based intervals.
///         Can be used in combination with `region`.
///
/// # Returns
///
/// A Polars dataframe.
///
/// # Example output
///
/// # Errors
/// If building of option-related structs fails, if BAM input
/// cannot be obtained, if preparing records fails, or running the
/// `nanalogue_core` commands fail
#[expect(
    clippy::too_many_arguments,
    clippy::fn_params_excessive_bools,
    reason = "python functions have many more args than rust, o.k. as some are optional"
)]
#[pyfunction]
#[pyo3(signature = (bam_path, treat_as_url = false, min_seq_len = 0, min_align_len = 0,
                    read_ids = HashSet::<String>::new(), threads = 2, include_zero_len = false,
                    read_filter = "", sample_fraction = 1.0, mapq_filter = 0,
                    exclude_mapq_unavail = false, region = "", full_region = false,
                    mod_strand = "", min_mod_qual = 0,
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

    // return output
    Ok(PyDataFrame(
        curr_reads_to_dataframe(df_collection.as_slice()).map_err(|e: Error| py_exception!(e))?,
    ))
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

    Ok(())
}
