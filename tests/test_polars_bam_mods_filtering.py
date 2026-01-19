# Testing filtering parameters for pynanalogue functions
# Uses dataclass pattern for managing many-parameter test cases

from dataclasses import dataclass, field, replace
from typing import Set
import pytest
import polars as pl
import pynanalogue


@dataclass
class TestInputOptions:
    """Test data builder for pynanalogue function parameters.

    Provides sensible defaults for all parameters, allowing tests to only
    override the specific parameters they need to test. This avoids the
    combinatorial explosion of testing all 19 parameters exhaustively.
    """

    bam_path: str = "test.bam"
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
    region: str = ""
    full_region: bool = False
    tag: str = ""
    mod_strand: str = ""
    min_mod_qual: int = 0
    reject_mod_qual_non_inclusive: tuple = (0, 0)
    trim_read_ends_mod: int = 0
    base_qual_filter_mod: int = 0
    mod_region: str = ""

    def as_dict(self):
        """Convert to dictionary for **kwargs unpacking"""
        return {k: v for k, v in self.__dict__.items()}


class TestInputBamFiltering:
    """Test filtering parameters related to InputBam struct"""

    def test_min_seq_len_filter(self, simple_bam):
        """Test that min_seq_len correctly filters reads"""
        base = TestInputOptions(bam_path=str(simple_bam))

        # Get all reads
        result_all = pynanalogue.polars_bam_mods(**base.as_dict())

        # Filter with min_seq_len=6000 (our test reads are 5000bp)
        params_filtered = replace(base, min_seq_len=6000)
        result_filtered = pynanalogue.polars_bam_mods(**params_filtered.as_dict())

        # Should filter out all reads since they're all 5kb
        assert len(result_all) > 0, "Expected some reads in unfiltered data"
        assert len(result_filtered) == 0, "Expected no reads with min_seq_len=6000"

    def test_min_align_len_filter(self, simple_bam):
        """Test that min_align_len correctly filters reads"""
        base = TestInputOptions(bam_path=str(simple_bam))

        # Get all reads
        result_all = pynanalogue.polars_bam_mods(**base.as_dict())

        # Filter with min_align_len=1 first and then with min_align_len=6000 (our test reads are 5000bp)
        params_filtered_1 = replace(base, min_align_len=1)
        result_filtered_1 = pynanalogue.polars_bam_mods(**params_filtered_1.as_dict())
        params_filtered_2 = replace(base, min_align_len=6000)
        result_filtered_2 = pynanalogue.polars_bam_mods(**params_filtered_2.as_dict())

        # Should filter out some reads as some are unmapped
        assert len(result_all) > 0, "Expected some reads in unfiltered data"
        assert len(result_filtered_1) < len(result_all), (
            "Unmapped reads must have filtered out"
        )

        # Should filter out all reads since they're all 5kb
        assert len(result_all) > 0, "Expected some reads in unfiltered data"
        assert len(result_filtered_2) == 0, "Expected no reads with min_align_len=6000"

    def test_mapq_filter_and_exclude_mapq_unavail(self, simple_bam):
        """Test that mapq_filter correctly filters reads"""
        base = TestInputOptions(bam_path=str(simple_bam))

        result_all = pynanalogue.polars_bam_mods(**base.as_dict())

        # Filter with very high mapq (our test data has mapq 10-20)
        # Note: Unmapped reads don't have MAPQ, so they pass through
        params_filtered = replace(base, mapq_filter=100)
        result_filtered = pynanalogue.polars_bam_mods(**params_filtered.as_dict())

        assert len(result_all) > 0
        assert len(result_filtered) < len(result_all)
        assert len(result_filtered) > 0  # Unmapped reads still present

        # Now exclude reads without mapq and verify we get zero results
        params_filtered_2 = replace(base, mapq_filter=100, exclude_mapq_unavail=True)
        result_filtered_2 = pynanalogue.polars_bam_mods(**params_filtered_2.as_dict())

        assert len(result_filtered_2) == 0, (
            "Expected no reads with mapq_filter=100 and exclude_mapq_unavail=True"
        )

    @pytest.mark.parametrize(
        "sample_fraction,expected_fraction",
        [
            (1.0, 1.0),  # No sampling
            (0.5, 0.5),  # Half sampling (approximate)
            (0.1, 0.1),  # 10% sampling (approximate)
        ],
    )
    def test_sample_fraction(self, simple_bam, sample_fraction, expected_fraction):
        """Test that sample_fraction approximately samples the expected proportion"""
        base = TestInputOptions(bam_path=str(simple_bam))

        # Get baseline count
        result_all = pynanalogue.polars_bam_mods(**base.as_dict())
        all_count = len(result_all.unique(subset=["read_id"]))

        # Sample
        params_sampled = replace(base, sample_fraction=sample_fraction)
        result_sampled = pynanalogue.polars_bam_mods(**params_sampled.as_dict())
        sampled_count = len(result_sampled.unique(subset=["read_id"]))

        if sample_fraction == 1.0:
            assert sampled_count == all_count
        else:
            # Allow 30% variance due to stochastic sampling
            expected = all_count * expected_fraction
            assert 0.7 * expected <= sampled_count <= 1.3 * expected, (
                f"Expected ~{expected} reads, got {sampled_count}"
            )

    def test_different_region_filters(self, simple_bam):
        """Test that region filtering works"""
        base = TestInputOptions(bam_path=str(simple_bam))

        # Test with a specific region (simulated BAM contigs are named contig_00000, contig_00001, etc.)
        params_region = replace(base, region="contig_00000")
        result = pynanalogue.polars_bam_mods(**params_region.as_dict())

        # Verify all results are from contig_00000
        if len(result) > 0:
            unique_contigs = result["contig"].unique().to_list()
            assert unique_contigs == ["contig_00000"], (
                f"Expected only contig_00000, got {unique_contigs}"
            )
        else:
            raise AssertionError("Expected some reads that map to contig_00000")

        # Test with the same region that if we request full region, then we don't get any reads
        # as no reads pass through the entire region.
        params_region_2 = replace(base, region="contig_00000", full_region=True)
        result_2 = pynanalogue.polars_bam_mods(**params_region_2.as_dict())
        assert len(result_2) == 0, (
            "No reads are expected to pass through the whole region"
        )

    def test_read_filter_primary_only(self, simple_bam):
        """Test that read_filter correctly filters by alignment type"""
        base = TestInputOptions(bam_path=str(simple_bam))

        # Filter to primary alignments only
        params_primary = replace(base, read_filter="primary_forward,primary_reverse")
        result_primary = pynanalogue.polars_bam_mods(**params_primary.as_dict())

        # Check that we only have primary alignments
        if len(result_primary) > 0:
            alignment_types = result_primary["alignment_type"].unique().to_list()
            assert all("primary" in atype for atype in alignment_types), (
                f"Expected only primary alignments, got {alignment_types}"
            )
        else:
            raise AssertionError("expected some primary reads!")

        # Test that whitespace-separated filter strings are trimmed and produce same results
        params_primary_2 = replace(base, read_filter="primary_forward, primary_reverse")
        result_primary_2 = pynanalogue.polars_bam_mods(**params_primary_2.as_dict())

        assert len(result_primary) == len(result_primary_2), (
            "Whitespace in comma-separated filter should be trimmed and produce same results"
        )

    def test_read_ids_filter_single(self, simple_bam):
        """Test filtering by a single read_id"""
        base = TestInputOptions(bam_path=str(simple_bam))

        # Load all data
        result_all = pynanalogue.polars_bam_mods(**base.as_dict())

        # Pick a random read id
        all_read_ids = result_all["read_id"].unique().to_list()
        selected_read_id = all_read_ids[0]

        # Count expected rows for this read id
        expected_count = len(result_all.filter(pl.col("read_id") == selected_read_id))

        # Filter by this single read id
        params_filtered = replace(base, read_ids={selected_read_id})
        result_filtered = pynanalogue.polars_bam_mods(**params_filtered.as_dict())

        assert len(result_filtered) == expected_count, (
            f"Expected {expected_count} rows for read_id {selected_read_id}, got {len(result_filtered)}"
        )
        assert result_filtered["read_id"].unique().to_list() == [selected_read_id]

    def test_read_ids_filter_two(self, simple_bam):
        """Test filtering by two read_ids"""
        base = TestInputOptions(bam_path=str(simple_bam))

        # Load all data
        result_all = pynanalogue.polars_bam_mods(**base.as_dict())

        # Pick two read ids
        all_read_ids = result_all["read_id"].unique().to_list()
        selected_read_ids = set(all_read_ids[:2])

        # Count expected rows for these read ids
        expected_count = len(
            result_all.filter(pl.col("read_id").is_in(selected_read_ids))
        )

        # Filter by these two read ids
        params_filtered = replace(base, read_ids=selected_read_ids)
        result_filtered = pynanalogue.polars_bam_mods(**params_filtered.as_dict())

        assert len(result_filtered) == expected_count, (
            f"Expected {expected_count} rows for read_ids {selected_read_ids}, got {len(result_filtered)}"
        )
        assert set(result_filtered["read_id"].unique().to_list()) == selected_read_ids


class TestInputModsFiltering:
    """Test filtering parameters related to InputMods struct"""

    def test_mod_strand_filter(self, simple_bam):
        """Test that mod_strand filtering works"""
        base = TestInputOptions(bam_path=str(simple_bam))

        result_all = pynanalogue.polars_bam_mods(**base.as_dict())

        # Filter to only basecalled complement strand
        params_bc_comp = replace(base, mod_strand="bc_comp")
        result_bc_comp = pynanalogue.polars_bam_mods(**params_bc_comp.as_dict())

        # Our test data has mods on basecalled strand, not complement
        assert len(result_all) > 0
        assert len(result_bc_comp) == 0, "Expected no mods on complement strand"

    def test_min_mod_qual_filter(self, simple_bam):
        """Test that min_mod_qual correctly filters low-quality mod calls"""
        base = TestInputOptions(bam_path=str(simple_bam))

        result_all = pynanalogue.polars_bam_mods(**base.as_dict())

        # Filter with high quality threshold
        params_high_qual = replace(base, min_mod_qual=200)
        result_high_qual = pynanalogue.polars_bam_mods(**params_high_qual.as_dict())

        # Should have fewer mods with higher quality threshold
        assert len(result_high_qual) < len(result_all)

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
        base = TestInputOptions(bam_path=str(simple_bam))

        params = replace(base, reject_mod_qual_non_inclusive=(low, high))

        if should_succeed:
            result = pynanalogue.polars_bam_mods(**params.as_dict())
            assert isinstance(result, pl.DataFrame)  # Should return DataFrame
        else:
            with pytest.raises(ValueError, match="low < high"):
                pynanalogue.polars_bam_mods(**params.as_dict())

    def test_trim_read_ends_mod(self, simple_bam):
        """Test that trim_read_ends_mod removes mods near read ends"""
        base = TestInputOptions(bam_path=str(simple_bam))

        result_all = pynanalogue.polars_bam_mods(**base.as_dict())

        # Trim 1000bp from each end
        params_trimmed = replace(base, trim_read_ends_mod=1000)
        result_trimmed = pynanalogue.polars_bam_mods(**params_trimmed.as_dict())

        # Should have fewer mods after trimming ends
        assert len(result_trimmed) < len(result_all)

    def test_base_qual_filter_mod(self, simple_bam):
        """Test that base_qual_filter_mod removes mods on low-quality bases"""
        base = TestInputOptions(bam_path=str(simple_bam))

        result_all = pynanalogue.polars_bam_mods(**base.as_dict())

        # Filter mods on bases with quality < 15
        params_qual_filtered = replace(base, base_qual_filter_mod=15)
        result_qual_filtered = pynanalogue.polars_bam_mods(
            **params_qual_filtered.as_dict()
        )

        # Should have fewer mods with quality filtering
        assert len(result_qual_filtered) < len(result_all)

    def test_mod_region_filter(self, simple_bam):
        """Test that mod_region filtering restricts results to specified regions"""
        base = TestInputOptions(bam_path=str(simple_bam))

        # Get all mods (no region filter)
        result_all = pynanalogue.polars_bam_mods(**base.as_dict())
        count_all = len(result_all)

        # Filter to just contig_00000
        params_contig = replace(base, mod_region="contig_00000")
        result_contig = pynanalogue.polars_bam_mods(**params_contig.as_dict())
        count_contig = len(result_contig)

        # Filter to specific range within contig_00000
        params_range = replace(base, mod_region="contig_00000:1000-2000")
        result_range = pynanalogue.polars_bam_mods(**params_range.as_dict())
        count_range = len(result_range)

        # Filtering by contig should give fewer results than no filter
        assert count_contig < count_all, (
            f"Expected fewer mods with contig filter, got {count_contig} vs {count_all}"
        )

        # Filtering by specific range should give even fewer results
        assert count_range < count_contig, (
            f"Expected fewer mods with range filter, got {count_range} vs {count_contig}"
        )

    def test_tag_filter(self, two_mods_bam):
        """Test that tag parameter filters to specific modification types.

        The two_mods_bam fixture has:
        - T modifications on minus strand (mod_code="T")
        - C modifications on plus strand (mod_code="76792")
        """
        base = TestInputOptions(bam_path=str(two_mods_bam))

        # Get all mods (no tag filter)
        result_all = pynanalogue.polars_bam_mods(**base.as_dict())
        mod_codes_all = set(result_all["mod_code"].unique().to_list())

        # Verify both mod types are present
        assert "T" in mod_codes_all, "Expected T mod in unfiltered data"
        assert "76792" in mod_codes_all, "Expected 76792 mod in unfiltered data"

        # Filter to only 76792 mods
        params_76792 = replace(base, tag="76792")
        result_76792 = pynanalogue.polars_bam_mods(**params_76792.as_dict())
        mod_codes_76792 = set(result_76792["mod_code"].unique().to_list())

        assert len(result_76792) > 0, "Expected some 76792 mods"
        assert mod_codes_76792 == {"76792"}, (
            f"Expected only 76792 mod, got {mod_codes_76792}"
        )

        # Filter to only T mods
        params_t = replace(base, tag="T")
        result_t = pynanalogue.polars_bam_mods(**params_t.as_dict())
        mod_codes_t = set(result_t["mod_code"].unique().to_list())

        assert len(result_t) > 0, "Expected some T mods"
        assert mod_codes_t == {"T"}, f"Expected only T mod, got {mod_codes_t}"

        # Verify counts add up
        assert len(result_76792) + len(result_t) == len(result_all), (
            f"Sum of filtered results ({len(result_76792)} + {len(result_t)}) "
            f"should equal unfiltered ({len(result_all)})"
        )
