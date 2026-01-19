# Testing filtering parameters for read_info function
# Uses dataclass pattern for managing many-parameter test cases

from dataclasses import dataclass, field, replace
from typing import Set
import json
import re
import pytest
import pynanalogue


def decode_read_info(result_bytes: bytes) -> list:
    """Decode the bytes output from read_info to a list of dicts"""
    return json.loads(result_bytes)


def parse_mod_count(mod_count_str: str) -> int:
    """Parse the mod_count string and extract the total count for the T+T mod.

    The format is like: "T+T:2104;(probabilities >= 0.5020, PHRED base qual >= 0)"
    or with multiple mods: "G-7200:0;T+T:3;(probabilities >= ...)"

    For our simple_bam test data, we only have T+T mods.
    """
    if mod_count_str == "NA":
        return 0

    # Match pattern like "T+T:2104" - looking for T+T specifically
    match = re.search(r"T\+T:(\d+)", mod_count_str)
    if match:
        return int(match.group(1))
    return 0


def get_total_mod_count(records: list) -> int:
    """Sum up all mod counts across all records"""
    total = 0
    for record in records:
        mod_count_str = record.get("mod_count", "NA")
        total += parse_mod_count(mod_count_str)
    return total


@dataclass
class InputOptions:
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
        base = InputOptions(bam_path=str(simple_bam))

        # Get all reads
        result_all = decode_read_info(pynanalogue.read_info(**base.as_dict()))

        # Filter with min_seq_len=6000 (our test reads are 5000bp)
        params_filtered = replace(base, min_seq_len=6000)
        result_filtered = decode_read_info(
            pynanalogue.read_info(**params_filtered.as_dict())
        )

        # Should filter out all reads since they're all 5kb
        assert len(result_all) > 0, "Expected some reads in unfiltered data"
        assert len(result_filtered) == 0, "Expected no reads with min_seq_len=6000"

    def test_min_align_len_filter(self, simple_bam):
        """Test that min_align_len correctly filters reads"""
        base = InputOptions(bam_path=str(simple_bam))

        # Get all reads
        result_all = decode_read_info(pynanalogue.read_info(**base.as_dict()))

        # Filter with min_align_len=1 first and then with min_align_len=6000
        params_filtered_1 = replace(base, min_align_len=1)
        result_filtered_1 = decode_read_info(
            pynanalogue.read_info(**params_filtered_1.as_dict())
        )
        params_filtered_2 = replace(base, min_align_len=6000)
        result_filtered_2 = decode_read_info(
            pynanalogue.read_info(**params_filtered_2.as_dict())
        )

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
        base = InputOptions(bam_path=str(simple_bam))

        result_all = decode_read_info(pynanalogue.read_info(**base.as_dict()))

        # Filter with very high mapq (our test data has mapq 10-20)
        # Note: Unmapped reads don't have MAPQ, so they pass through
        params_filtered = replace(base, mapq_filter=100)
        result_filtered = decode_read_info(
            pynanalogue.read_info(**params_filtered.as_dict())
        )

        assert len(result_all) > 0
        assert len(result_filtered) < len(result_all)
        assert len(result_filtered) > 0  # Unmapped reads still present

        # Now exclude reads without mapq and verify we get zero results
        params_filtered_2 = replace(base, mapq_filter=100, exclude_mapq_unavail=True)
        result_filtered_2 = decode_read_info(
            pynanalogue.read_info(**params_filtered_2.as_dict())
        )

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
        base = InputOptions(bam_path=str(simple_bam))

        # Get baseline count
        result_all = decode_read_info(pynanalogue.read_info(**base.as_dict()))
        all_read_ids = {r["read_id"] for r in result_all}
        all_count = len(all_read_ids)

        # Sample
        params_sampled = replace(base, sample_fraction=sample_fraction)
        result_sampled = decode_read_info(
            pynanalogue.read_info(**params_sampled.as_dict())
        )
        sampled_read_ids = {r["read_id"] for r in result_sampled}
        sampled_count = len(sampled_read_ids)

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
        base = InputOptions(bam_path=str(simple_bam))

        # Test with a specific region (simulated BAM contigs are named contig_00000, etc.)
        params_region = replace(base, region="contig_00000")
        result = decode_read_info(pynanalogue.read_info(**params_region.as_dict()))

        # Verify all results are from contig_00000
        if len(result) > 0:
            unique_contigs = {r.get("contig") for r in result if r.get("contig")}
            assert unique_contigs == {"contig_00000"}, (
                f"Expected only contig_00000, got {unique_contigs}"
            )
        else:
            raise AssertionError("Expected some reads that map to contig_00000")

        # Test with the same region that if we request full region, then we don't get
        # any reads as no reads pass through the entire region.
        params_region_2 = replace(base, region="contig_00000", full_region=True)
        result_2 = decode_read_info(pynanalogue.read_info(**params_region_2.as_dict()))
        assert len(result_2) == 0, (
            "No reads are expected to pass through the whole region"
        )

    def test_read_filter_primary_only(self, simple_bam):
        """Test that read_filter correctly filters by alignment type"""
        base = InputOptions(bam_path=str(simple_bam))

        # Filter to primary alignments only
        params_primary = replace(base, read_filter="primary_forward,primary_reverse")
        result_primary = decode_read_info(
            pynanalogue.read_info(**params_primary.as_dict())
        )

        # Check that we only have primary alignments
        if len(result_primary) > 0:
            alignment_types = {r["alignment_type"] for r in result_primary}
            assert all("primary" in atype for atype in alignment_types), (
                f"Expected only primary alignments, got {alignment_types}"
            )
        else:
            raise AssertionError("expected some primary reads!")

        # Test that whitespace-separated filter strings are trimmed and produce same results
        params_primary_2 = replace(base, read_filter="primary_forward, primary_reverse")
        result_primary_2 = decode_read_info(
            pynanalogue.read_info(**params_primary_2.as_dict())
        )

        assert len(result_primary) == len(result_primary_2), (
            "Whitespace in comma-separated filter should be trimmed and produce same results"
        )

    def test_read_ids_filter_single(self, simple_bam):
        """Test filtering by a single read_id"""
        base = InputOptions(bam_path=str(simple_bam))

        # Load all data
        result_all = decode_read_info(pynanalogue.read_info(**base.as_dict()))

        # Pick a random read id
        all_read_ids = list({r["read_id"] for r in result_all})
        selected_read_id = all_read_ids[0]

        # Count expected records for this read id
        expected_count = sum(1 for r in result_all if r["read_id"] == selected_read_id)

        # Filter by this single read id
        params_filtered = replace(base, read_ids={selected_read_id})
        result_filtered = decode_read_info(
            pynanalogue.read_info(**params_filtered.as_dict())
        )

        assert len(result_filtered) == expected_count, (
            f"Expected {expected_count} records for read_id {selected_read_id}, "
            f"got {len(result_filtered)}"
        )
        assert all(r["read_id"] == selected_read_id for r in result_filtered)

    def test_read_ids_filter_two(self, simple_bam):
        """Test filtering by two read_ids"""
        base = InputOptions(bam_path=str(simple_bam))

        # Load all data
        result_all = decode_read_info(pynanalogue.read_info(**base.as_dict()))

        # Pick two read ids
        all_read_ids = list({r["read_id"] for r in result_all})
        selected_read_ids = set(all_read_ids[:2])

        # Count expected records for these read ids
        expected_count = sum(1 for r in result_all if r["read_id"] in selected_read_ids)

        # Filter by these two read ids
        params_filtered = replace(base, read_ids=selected_read_ids)
        result_filtered = decode_read_info(
            pynanalogue.read_info(**params_filtered.as_dict())
        )

        assert len(result_filtered) == expected_count, (
            f"Expected {expected_count} records for read_ids {selected_read_ids}, "
            f"got {len(result_filtered)}"
        )
        assert {r["read_id"] for r in result_filtered} == selected_read_ids


class TestInputModsFiltering:
    """Test filtering parameters related to InputMods struct"""

    def test_mod_strand_filter(self, simple_bam):
        """Test that mod_strand filtering works"""
        base = InputOptions(bam_path=str(simple_bam))

        result_all = decode_read_info(pynanalogue.read_info(**base.as_dict()))
        total_mods_all = get_total_mod_count(result_all)

        # Filter to only basecalled complement strand
        params_bc_comp = replace(base, mod_strand="bc_comp")
        result_bc_comp = decode_read_info(
            pynanalogue.read_info(**params_bc_comp.as_dict())
        )
        total_mods_bc_comp = get_total_mod_count(result_bc_comp)

        # Our test data has mods on basecalled strand, not complement
        # Same number of records but zero mod counts on complement
        assert len(result_all) == len(result_bc_comp), "Expected same number of records"
        assert total_mods_all > 0, "Expected some mods in unfiltered data"
        assert total_mods_bc_comp == 0, "Expected no mods on complement strand"

    def test_min_mod_qual_filter(self, simple_bam):
        """Test that min_mod_qual correctly filters low-quality mod calls"""
        base = InputOptions(bam_path=str(simple_bam))

        result_all = decode_read_info(pynanalogue.read_info(**base.as_dict()))
        total_mods_all = get_total_mod_count(result_all)

        # Filter with high quality threshold
        params_high_qual = replace(base, min_mod_qual=200)
        result_high_qual = decode_read_info(
            pynanalogue.read_info(**params_high_qual.as_dict())
        )
        total_mods_high_qual = get_total_mod_count(result_high_qual)

        # Same number of records but fewer mods with higher quality threshold
        assert len(result_all) == len(result_high_qual), (
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
        base = InputOptions(bam_path=str(simple_bam))

        params = replace(base, reject_mod_qual_non_inclusive=(low, high))

        if should_succeed:
            result = pynanalogue.read_info(**params.as_dict())
            decoded = decode_read_info(result)
            assert isinstance(decoded, list)  # Should return valid JSON list
        else:
            with pytest.raises(ValueError, match="low < high"):
                pynanalogue.read_info(**params.as_dict())

    def test_trim_read_ends_mod(self, simple_bam):
        """Test that trim_read_ends_mod removes mods near read ends"""
        base = InputOptions(bam_path=str(simple_bam))

        result_all = decode_read_info(pynanalogue.read_info(**base.as_dict()))
        total_mods_all = get_total_mod_count(result_all)

        # Trim 1000bp from each end
        params_trimmed = replace(base, trim_read_ends_mod=1000)
        result_trimmed = decode_read_info(
            pynanalogue.read_info(**params_trimmed.as_dict())
        )
        total_mods_trimmed = get_total_mod_count(result_trimmed)

        # Same number of records but fewer mods after trimming ends
        assert len(result_all) == len(result_trimmed), "Expected same number of records"
        assert total_mods_trimmed < total_mods_all, (
            f"Expected fewer mods after trimming, "
            f"got {total_mods_trimmed} vs {total_mods_all}"
        )

    def test_base_qual_filter_mod(self, simple_bam):
        """Test that base_qual_filter_mod removes mods on low-quality bases"""
        base = InputOptions(bam_path=str(simple_bam))

        result_all = decode_read_info(pynanalogue.read_info(**base.as_dict()))
        total_mods_all = get_total_mod_count(result_all)

        # Filter mods on bases with quality < 15
        params_qual_filtered = replace(base, base_qual_filter_mod=15)
        result_qual_filtered = decode_read_info(
            pynanalogue.read_info(**params_qual_filtered.as_dict())
        )
        total_mods_qual_filtered = get_total_mod_count(result_qual_filtered)

        # Same number of records but fewer mods with quality filtering
        assert len(result_all) == len(result_qual_filtered), (
            "Expected same number of records"
        )
        assert total_mods_qual_filtered < total_mods_all, (
            f"Expected fewer mods with quality filter, "
            f"got {total_mods_qual_filtered} vs {total_mods_all}"
        )

    def test_mod_region_filter(self, simple_bam):
        """Test that mod_region filtering restricts results to specified regions"""
        base = InputOptions(bam_path=str(simple_bam))

        # Get all mods (no region filter)
        result_all = decode_read_info(pynanalogue.read_info(**base.as_dict()))
        total_mods_all = get_total_mod_count(result_all)

        # Filter to just contig_00000
        params_contig = replace(base, mod_region="contig_00000")
        result_contig = decode_read_info(
            pynanalogue.read_info(**params_contig.as_dict())
        )
        total_mods_contig = get_total_mod_count(result_contig)

        # Filter to specific range within contig_00000
        params_range = replace(base, mod_region="contig_00000:1000-2000")
        result_range = decode_read_info(pynanalogue.read_info(**params_range.as_dict()))
        total_mods_range = get_total_mod_count(result_range)

        # Same number of records in all cases
        assert len(result_all) == len(result_contig) == len(result_range), (
            "Expected same number of records"
        )

        # Filtering by contig should give fewer mods than no filter
        assert total_mods_contig < total_mods_all, (
            f"Expected fewer mods with contig filter, "
            f"got {total_mods_contig} vs {total_mods_all}"
        )

        # Filtering by specific range should give even fewer mods
        assert total_mods_range < total_mods_contig, (
            f"Expected fewer mods with range filter, "
            f"got {total_mods_range} vs {total_mods_contig}"
        )
