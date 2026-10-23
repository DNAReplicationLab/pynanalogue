# Tests for pynanalogue.seq_table() function
# Validates sequence extraction from BAM files with modification markers

from pathlib import Path

import polars as pl
import pynanalogue
import pytest

# Test data paths
TEST_DIR = Path(__file__).parent / "data" / "examples"
BAM_FORWARD = str(TEST_DIR / "example_pynanalogue_1.bam")
BAM_FIRST_READ_REV = str(TEST_DIR / "example_pynanalogue_1_first_read_rev.bam")


@pytest.mark.parametrize(
    "region,expected_rows",
    [
        pytest.param("dummyI:9-17", 1, id="exact-alignment"),
        pytest.param("dummyI:10-16", 1, id="inside-alignment"),
        pytest.param("dummyI:8-17", 0, id="extends-left"),
        pytest.param("dummyI:9-18", 0, id="extends-right"),
        pytest.param("dummyI:9-", 0, id="open-ended"),
        pytest.param("dummyI", 0, id="whole-contig"),
    ],
)
def test_seq_table_requires_full_region_and_preserves_schema(region, expected_rows):
    """Only fully spanning alignments are retained, with a stable empty schema."""
    # This fixture's only dummyI alignment spans [9, 17).
    bam = str(TEST_DIR / "example_1.bam")
    df = pynanalogue.seq_table(bam, region=region)
    assert df.height == expected_rows
    assert df["read_id"].to_list() == (
        ["5d10eb9a-aae1-4db8-8ec6-7ebb34d32575"] if expected_rows else []
    )
    assert df.columns == ["read_id", "sequence", "qualities"]
    assert df.schema == {
        "read_id": pl.String,
        "sequence": pl.String,
        "qualities": pl.String,
    }


@pytest.mark.parametrize(
    "bam_file,expected_0_seq,expected_1_seq",
    [
        (BAM_FORWARD, "AZGTAZGTAZ", "ACGTACGTAC"),
        (BAM_FIRST_READ_REV, "ACZTACZTAC", "ACGTACGTAC"),
    ],
    ids=["forward", "first_read_rev"],
)
def test_seq_table_region_0_10(bam_file, expected_0_seq, expected_1_seq):
    """Test seq_table on region contig_00000:0-10"""
    df = pynanalogue.seq_table(bam_path=bam_file, region="contig_00000:0-10")

    assert isinstance(df, pl.DataFrame)
    assert len(df) == 2

    expected = {
        "0.": {
            "sequence": expected_0_seq,
            "qualities": "20.20.20.20.20.20.20.20.20.20",
        },
        "1.": {
            "sequence": expected_1_seq,
            "qualities": "30.30.30.30.30.30.30.30.30.30",
        },
    }

    for row in df.iter_rows(named=True):
        read_id = row["read_id"]
        prefix = read_id[:2]
        assert prefix in expected, f"Unexpected read_id prefix: {prefix}"
        assert row["sequence"] == expected[prefix]["sequence"], (
            f"Sequence mismatch for {prefix}"
        )
        assert row["qualities"] == expected[prefix]["qualities"], (
            f"Qualities mismatch for {prefix}"
        )


@pytest.mark.parametrize(
    "bam_file,expected_0_seq,expected_1_seq",
    [
        (BAM_FORWARD, "TAZGT.....", "TACGT....."),
        (BAM_FIRST_READ_REV, "TACZT.....", "TACGT....."),
    ],
    ids=["forward", "first_read_rev"],
)
def test_seq_table_region_15_25(bam_file, expected_0_seq, expected_1_seq):
    """Test seq_table on region contig_00000:15-25 with partial coverage"""
    df = pynanalogue.seq_table(bam_path=bam_file, region="contig_00000:15-25")

    assert isinstance(df, pl.DataFrame)
    assert len(df) == 2

    expected = {
        "0.": {
            "sequence": expected_0_seq,
            "qualities": "20.20.20.20.20.255.255.255.255.255",
        },
        "1.": {
            "sequence": expected_1_seq,
            "qualities": "30.30.30.30.30.255.255.255.255.255",
        },
    }

    for row in df.iter_rows(named=True):
        read_id = row["read_id"]
        prefix = read_id[:2]
        assert prefix in expected, f"Unexpected read_id prefix: {prefix}"
        assert row["sequence"] == expected[prefix]["sequence"], (
            f"Sequence mismatch for {prefix}"
        )
        assert row["qualities"] == expected[prefix]["qualities"], (
            f"Qualities mismatch for {prefix}"
        )


@pytest.mark.parametrize(
    "bam_file,expected_0_seq,expected_1_seq",
    [
        (BAM_FORWARD, "TAZGTztztAZGTA", "TACGTctctACGTA"),
        (BAM_FIRST_READ_REV, "TACZTctctACZTA", "TACGTctctACGTA"),
    ],
    ids=["forward", "first_read_rev"],
)
def test_seq_table_region_95_105(bam_file, expected_0_seq, expected_1_seq):
    """Test seq_table on region contig_00000:95-105 with insertions"""
    df = pynanalogue.seq_table(bam_path=bam_file, region="contig_00000:95-105")

    assert isinstance(df, pl.DataFrame)
    assert len(df) == 2

    expected = {
        "0.": {
            "sequence": expected_0_seq,
            "qualities": "20.20.20.20.20.20.20.20.20.20.20.20.20.20",
        },
        "1.": {
            "sequence": expected_1_seq,
            "qualities": "30.30.30.30.30.30.30.30.30.30.30.30.30.30",
        },
    }

    for row in df.iter_rows(named=True):
        read_id = row["read_id"]
        prefix = read_id[:2]
        assert prefix in expected, f"Unexpected read_id prefix: {prefix}"
        assert row["sequence"] == expected[prefix]["sequence"], (
            f"Sequence mismatch for {prefix}"
        )
        assert row["qualities"] == expected[prefix]["qualities"], (
            f"Qualities mismatch for {prefix}"
        )


@pytest.mark.parametrize(
    "bam_file,expected_0_seq,expected_1_seq",
    [
        (BAM_FORWARD, "GTAZGTAZGT", "GTACGTACGT"),
        (BAM_FIRST_READ_REV, "ZTACZTACZT", "GTACGTACGT"),
    ],
    ids=["forward", "first_read_rev"],
)
def test_seq_table_region_190_200(bam_file, expected_0_seq, expected_1_seq):
    """Test seq_table on region contig_00000:190-200"""
    df = pynanalogue.seq_table(bam_path=bam_file, region="contig_00000:190-200")

    assert isinstance(df, pl.DataFrame)
    assert len(df) == 2

    expected = {
        "0.": {
            "sequence": expected_0_seq,
            "qualities": "20.20.20.20.20.20.20.20.20.20",
        },
        "1.": {
            "sequence": expected_1_seq,
            "qualities": "30.30.30.30.30.30.30.30.30.30",
        },
    }

    for row in df.iter_rows(named=True):
        read_id = row["read_id"]
        prefix = read_id[:2]
        assert prefix in expected, f"Unexpected read_id prefix: {prefix}"
        assert row["sequence"] == expected[prefix]["sequence"], (
            f"Sequence mismatch for {prefix}"
        )
        assert row["qualities"] == expected[prefix]["qualities"], (
            f"Qualities mismatch for {prefix}"
        )
