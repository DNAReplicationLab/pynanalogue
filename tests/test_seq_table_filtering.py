# Testing filtering parameters for seq_table function
# Uses dataclass pattern for managing many-parameter test cases

from dataclasses import dataclass, field, replace
from typing import Set
import json
import polars as pl
import pytest
import pynanalogue


def get_row_count(df: pl.DataFrame) -> int:
    """Get the number of rows in the DataFrame"""
    return len(df)


def count_mods_in_sequence(sequence: str) -> int:
    """Count modification markers (Z and z) in a sequence string.

    Args:
        sequence: Sequence string where Z = mod on reference, z = mod in insertion

    Returns:
        Total count of Z and z characters
    """
    return sequence.count("Z") + sequence.count("z")


def get_total_mod_count(df: pl.DataFrame) -> int:
    """Sum up all mod counts across all sequences in the DataFrame.

    Args:
        df: DataFrame with a 'sequence' column

    Returns:
        Total count of modifications (Z and z) across all rows
    """
    total = 0
    for seq in df["sequence"].to_list():
        total += count_mods_in_sequence(seq)
    return total


@dataclass
class SeqTableInputOptions:
    """Test data builder for pynanalogue.seq_table function parameters.

    Provides sensible defaults for all parameters, allowing tests to only
    override the specific parameters they need to test.

    Note: seq_table requires 'region' parameter and has full_region=True
    hardcoded internally. It also doesn't have mod_region parameter.
    """

    bam_path: str = "test.bam"
    region: str = "contig_00000:4000-6000"
    treat_as_url: bool = False
    min_seq_len: int = 0
    min_align_len: int = 0
    read_ids: Set[str] = field(default_factory=set)
    threads: int = 2
    include_zero_len: bool = False
    read_filter: str = ""
    sample_fraction: float = 1.0
    mapq_filter: int = 0
    exclude_mapq_unavail: bool = False
    tag: str = ""
    mod_strand: str = ""
    min_mod_qual: int = 0
    reject_mod_qual_non_inclusive: tuple = (0, 0)
    trim_read_ends_mod: int = 0
    base_qual_filter_mod: int = 0

    def as_dict(self):
        """Convert to dictionary for **kwargs unpacking"""
        return {k: v for k, v in self.__dict__.items()}


class TestInputBamFiltering:
    """Test filtering parameters related to InputBam struct for seq_table"""

    def test_min_seq_len_filter(self, simple_bam):
        """Test that min_seq_len correctly filters reads"""
        base = SeqTableInputOptions(bam_path=str(simple_bam))

        # Get all reads
        result_all = pynanalogue.seq_table(**base.as_dict())
        all_count = get_row_count(result_all)

        # Filter with min_seq_len=6000 (our test reads are 5000bp)
        params_filtered = replace(base, min_seq_len=6000)
        result_filtered = pynanalogue.seq_table(**params_filtered.as_dict())
        filtered_count = get_row_count(result_filtered)

        # Should filter out all reads since they're all 5kb
        assert all_count > 0, "Expected some reads in unfiltered data"
        assert filtered_count == 0, "Expected no reads with min_seq_len=6000"

    def test_min_align_len_filter(self, simple_bam):
        """Test that min_align_len correctly filters reads"""
        base = SeqTableInputOptions(bam_path=str(simple_bam))

        # Get all reads
        result_all = pynanalogue.seq_table(**base.as_dict())
        all_count = get_row_count(result_all)

        # Filter with min_align_len=6000 (reads are 5000bp)
        params_filtered = replace(base, min_align_len=6000)
        result_filtered = pynanalogue.seq_table(**params_filtered.as_dict())
        filtered_count = get_row_count(result_filtered)

        # Should filter out all reads since they're all 5kb
        assert all_count > 0, "Expected some reads in unfiltered data"
        assert filtered_count == 0, "Expected no reads with min_align_len=6000"

    def test_mapq_filter(self, simple_bam):
        """Test that mapq_filter correctly filters reads"""
        base = SeqTableInputOptions(bam_path=str(simple_bam))

        result_all = pynanalogue.seq_table(**base.as_dict())
        all_count = get_row_count(result_all)

        # Filter with very high mapq (our test data has mapq 10-20)
        # Note: seq_table has full_region=True, so unmapped reads are already excluded
        params_filtered = replace(base, mapq_filter=100)
        result_filtered = pynanalogue.seq_table(**params_filtered.as_dict())
        filtered_count = get_row_count(result_filtered)

        assert all_count > 0
        assert filtered_count == 0, "Expected no reads with mapq_filter=100"

    @pytest.mark.parametrize(
        "sample_fraction,expected_fraction",
        [
            (1.0, 1.0),  # No sampling
            (0.5, 0.5),  # Half sampling (approximate)
            # (0.1, 0.1), 10% sampling (approximate)
            # We are removing this test case due to fluctuations.
            # We are now querying all reads that pass fully through
            # a small region, so the number of reads we are dealing with
            # is small and a coin toss with 10% probability of success
            # fluctuates greatly when sample size is small.
        ],
    )
    def test_sample_fraction(self, simple_bam, sample_fraction, expected_fraction):
        """Test that sample_fraction approximately samples the expected proportion"""
        base = SeqTableInputOptions(bam_path=str(simple_bam))

        # Get baseline count
        result_all = pynanalogue.seq_table(**base.as_dict())
        all_count = get_row_count(result_all)

        # Sample
        params_sampled = replace(base, sample_fraction=sample_fraction)
        result_sampled = pynanalogue.seq_table(**params_sampled.as_dict())
        sampled_count = get_row_count(result_sampled)

        if sample_fraction == 1.0:
            assert sampled_count == all_count
        else:
            # Allow 30% variance due to stochastic sampling
            expected = all_count * expected_fraction
            assert 0.7 * expected <= sampled_count <= 1.3 * expected, (
                f"Expected ~{expected} reads, got {sampled_count}"
            )

    def test_read_filter_primary_only(self, simple_bam):
        """Test that read_filter correctly filters by alignment type"""
        base = SeqTableInputOptions(bam_path=str(simple_bam))

        # Get all reads
        result_all = pynanalogue.seq_table(**base.as_dict())
        all_count = get_row_count(result_all)

        # Filter to primary alignments only
        params_primary = replace(base, read_filter="primary_forward,primary_reverse")
        result_primary = pynanalogue.seq_table(**params_primary.as_dict())
        primary_count = get_row_count(result_primary)

        # Should have at least some primary reads
        assert primary_count > 0, "Expected some primary reads"
        # Primary count should be less than all (some may be secondary/supplementary)
        assert primary_count < all_count

        # Test that whitespace-separated filter strings are trimmed and produce same results
        params_primary_2 = replace(base, read_filter="primary_forward, primary_reverse")
        result_primary_2 = pynanalogue.seq_table(**params_primary_2.as_dict())
        primary_count_2 = get_row_count(result_primary_2)

        assert primary_count == primary_count_2, (
            "Whitespace in comma-separated filter should be trimmed and produce same results"
        )

    def test_read_ids_filter_single(self, simple_bam):
        """Test filtering by a single read_id"""
        base = SeqTableInputOptions(bam_path=str(simple_bam))

        # Load all data
        result_all = pynanalogue.seq_table(**base.as_dict())

        # Pick a random read id from the results
        all_read_ids = result_all["read_id"].unique().to_list()
        if len(all_read_ids) == 0:
            pytest.skip("No reads in result to test with")

        selected_read_id = all_read_ids[0]

        # Count expected records for this read id
        expected_count = result_all.filter(pl.col("read_id") == selected_read_id).height

        # Filter by this single read id
        params_filtered = replace(base, read_ids={selected_read_id})
        result_filtered = pynanalogue.seq_table(**params_filtered.as_dict())
        filtered_count = get_row_count(result_filtered)

        assert filtered_count == expected_count, (
            f"Expected {expected_count} records for read_id {selected_read_id}, "
            f"got {filtered_count}"
        )
        assert all(
            rid == selected_read_id for rid in result_filtered["read_id"].to_list()
        )

    def test_read_ids_filter_two(self, simple_bam):
        """Test filtering by two read_ids"""
        base = SeqTableInputOptions(bam_path=str(simple_bam))

        # Load all data
        result_all = pynanalogue.seq_table(**base.as_dict())

        # Pick two read ids
        all_read_ids = result_all["read_id"].unique().to_list()
        if len(all_read_ids) < 2:
            pytest.skip("Not enough reads to test with two read_ids")

        selected_read_ids = set(all_read_ids[:2])

        # Count expected records for these read ids
        expected_count = result_all.filter(
            pl.col("read_id").is_in(list(selected_read_ids))
        ).height

        # Filter by these two read ids
        params_filtered = replace(base, read_ids=selected_read_ids)
        result_filtered = pynanalogue.seq_table(**params_filtered.as_dict())
        filtered_count = get_row_count(result_filtered)

        assert filtered_count == expected_count, (
            f"Expected {expected_count} records for read_ids {selected_read_ids}, "
            f"got {filtered_count}"
        )
        assert set(result_filtered["read_id"].to_list()) == selected_read_ids


class TestInputModsFiltering:
    """Test filtering parameters related to InputMods struct for seq_table"""

    def test_mod_strand_filter(self, simple_bam):
        """Test that mod_strand filtering works"""
        base = SeqTableInputOptions(bam_path=str(simple_bam))

        result_all = pynanalogue.seq_table(**base.as_dict())
        total_mods_all = get_total_mod_count(result_all)

        # Filter to only basecalled complement strand
        params_bc_comp = replace(base, mod_strand="bc_comp")
        result_bc_comp = pynanalogue.seq_table(**params_bc_comp.as_dict())
        total_mods_bc_comp = get_total_mod_count(result_bc_comp)

        # Our test data has mods on basecalled strand, not complement
        # Same number of records but zero mod counts on complement
        assert get_row_count(result_all) == get_row_count(result_bc_comp), (
            "Expected same number of records"
        )
        assert total_mods_all > 0, "Expected some mods in unfiltered data"
        assert total_mods_bc_comp == 0, "Expected no mods on complement strand"

    def test_min_mod_qual_filter(self, simple_bam):
        """Test that min_mod_qual correctly filters low-quality mod calls"""
        base = SeqTableInputOptions(bam_path=str(simple_bam))

        result_all = pynanalogue.seq_table(**base.as_dict())
        total_mods_all = get_total_mod_count(result_all)

        # Filter with high quality threshold
        params_high_qual = replace(base, min_mod_qual=200)
        result_high_qual = pynanalogue.seq_table(**params_high_qual.as_dict())
        total_mods_high_qual = get_total_mod_count(result_high_qual)

        # Same number of records but fewer mods with higher quality threshold
        assert get_row_count(result_all) == get_row_count(result_high_qual), (
            "Expected same number of records"
        )
        assert total_mods_high_qual < total_mods_all, (
            f"Expected fewer mods with high quality filter, "
            f"got {total_mods_high_qual} vs {total_mods_all}"
        )

    @pytest.mark.parametrize(
        "low,high,should_succeed",
        [
            (0, 0, True),  # Equal: should succeed with GtEq variant
            (0, 1, True),  # Diff of 1: should succeed with GtEq variant
            (0, 3, True),  # Diff > 1: should succeed with Both variant
            (100, 200, True),  # Valid range
            (200, 100, False),  # Invalid: low > high, should fail
        ],
    )
    def test_reject_mod_qual_validation(self, simple_bam, low, high, should_succeed):
        """Test reject_mod_qual_non_inclusive validation logic"""
        base = SeqTableInputOptions(bam_path=str(simple_bam))

        params = replace(base, reject_mod_qual_non_inclusive=(low, high))

        if should_succeed:
            result = pynanalogue.seq_table(**params.as_dict())
            assert isinstance(result, pl.DataFrame)
        else:
            with pytest.raises(ValueError, match="low < high"):
                pynanalogue.seq_table(**params.as_dict())

    def test_trim_read_ends_mod(self, simple_bam):
        """Test that trim_read_ends_mod removes mods near read ends"""
        base = SeqTableInputOptions(bam_path=str(simple_bam))

        result_all = pynanalogue.seq_table(**base.as_dict())
        total_mods_all = get_total_mod_count(result_all)

        # Trim 1000bp from each end
        params_trimmed = replace(base, trim_read_ends_mod=1000)
        result_trimmed = pynanalogue.seq_table(**params_trimmed.as_dict())
        total_mods_trimmed = get_total_mod_count(result_trimmed)

        # Same number of records but fewer mods after trimming ends
        assert get_row_count(result_all) == get_row_count(result_trimmed), (
            "Expected same number of records"
        )
        assert total_mods_trimmed < total_mods_all, (
            f"Expected fewer mods after trimming, "
            f"got {total_mods_trimmed} vs {total_mods_all}"
        )

    def test_base_qual_filter_mod(self, simple_bam):
        """Test that base_qual_filter_mod removes mods on low-quality bases"""
        base = SeqTableInputOptions(bam_path=str(simple_bam))

        result_all = pynanalogue.seq_table(**base.as_dict())
        total_mods_all = get_total_mod_count(result_all)

        # Filter mods on bases with quality < 15
        params_qual_filtered = replace(base, base_qual_filter_mod=15)
        result_qual_filtered = pynanalogue.seq_table(**params_qual_filtered.as_dict())
        total_mods_qual_filtered = get_total_mod_count(result_qual_filtered)

        # Same number of records but fewer mods with quality filtering
        assert get_row_count(result_all) == get_row_count(result_qual_filtered), (
            "Expected same number of records"
        )
        assert total_mods_qual_filtered < total_mods_all, (
            f"Expected fewer mods with quality filter, "
            f"got {total_mods_qual_filtered} vs {total_mods_all}"
        )


class TestTagFiltering:
    """Test tag parameter filtering for seq_table"""

    def test_tag_filter(self, two_mods_bam):
        """Test that tag parameter filters to specific modification types.

        The two_mods_bam fixture has:
        - T modifications on minus strand
        - C modifications on plus strand (ChEBI code 76792)
        """
        base = SeqTableInputOptions(bam_path=str(two_mods_bam))

        # Get all mods (no tag filter)
        result_all = pynanalogue.seq_table(**base.as_dict())
        total_mods_all = get_total_mod_count(result_all)

        # Filter to only 76792 mods
        params_76792 = replace(base, tag="76792")
        result_76792 = pynanalogue.seq_table(**params_76792.as_dict())
        total_mods_76792 = get_total_mod_count(result_76792)

        # Filter to only T mods
        params_t = replace(base, tag="T")
        result_t = pynanalogue.seq_table(**params_t.as_dict())
        total_mods_t = get_total_mod_count(result_t)

        # Total mods when filtering by specific tag should be less than all mods
        assert total_mods_all > 0, "Expected some mods in unfiltered data"

        # The sum of mods from each tag filter should equal total mods
        # (since mods don't overlap in our test data)
        assert total_mods_76792 + total_mods_t == total_mods_all, (
            f"Expected 76792 mods ({total_mods_76792}) + T mods ({total_mods_t}) "
            f"to equal total mods ({total_mods_all})"
        )

        # Each tag filter should have some mods (both are present in test data)
        assert total_mods_76792 > 0, "Expected some 76792 mods"
        assert total_mods_t > 0, "Expected some T mods"


class TestExcludeMapqUnavail:
    """Test exclude_mapq_unavail flag with custom BAM having MAPQ=255"""

    def test_exclude_mapq_unavail(self, tmp_path):
        """Test that exclude_mapq_unavail correctly filters reads with MAPQ=255.

        MAPQ=255 indicates mapping quality is unavailable. This test creates a
        BAM with all reads having MAPQ=255, then verifies:
        - Without exclude_mapq_unavail, reads pass through (MAPQ 255 passes any filter)
        - With exclude_mapq_unavail=True, all reads are filtered out
        """
        # Create BAM with MAPQ=255 (unavailable)
        config = {
            "contigs": {"number": 1, "len_range": [10000, 10000]},
            "reads": [
                {
                    "number": 1000,
                    "mapq_range": [255, 255],
                    "base_qual_range": [20, 30],
                    "len_range": [0.5, 0.5],
                    "insert_middle": "ATCG",
                    "mods": [],
                }
            ],
        }

        bam_path = tmp_path / "mapq_unavail.bam"
        fasta_path = tmp_path / "mapq_unavail.fasta"
        pynanalogue.simulate_mod_bam(
            json_config=json.dumps(config),
            bam_path=str(bam_path),
            fasta_path=str(fasta_path),
        )

        base = SeqTableInputOptions(
            bam_path=str(bam_path), region="contig_00000:4000-6000"
        )

        # Without exclude_mapq_unavail, reads with MAPQ=255 should pass through
        result_without_flag = pynanalogue.seq_table(**base.as_dict())
        count_without_flag = get_row_count(result_without_flag)

        # With exclude_mapq_unavail=True, all reads should be filtered out
        params_with_flag = replace(base, exclude_mapq_unavail=True)
        result_with_flag = pynanalogue.seq_table(**params_with_flag.as_dict())
        count_with_flag = get_row_count(result_with_flag)

        assert count_without_flag > 0, (
            "Expected some reads without exclude_mapq_unavail flag"
        )
        assert count_with_flag == 0, "Expected no reads with exclude_mapq_unavail=True"
