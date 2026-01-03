# Tests for pynanalogue.polars_bam_mods() function
# Validates Polars DataFrame output, filtering, and integration with example data

import pytest
from pathlib import Path
import polars as pl
import pynanalogue


def _test_polars_bam_mods_against_expected(example_number):
    """Helper function to test polars_bam_mods output against expected TSV output for a given example

    Args:
        example_number: The example number (e.g., 1, 3, 7)
    """
    # Paths to test data
    test_dir = Path(__file__).parent
    bam_path = test_dir / "data" / "examples" / f"example_{example_number}.bam"
    expected_path = (
        test_dir
        / "data"
        / "expected_outputs"
        / f"example_{example_number}_polars_bam_mods.tsv"
    )

    # Run polars_bam_mods
    result_df = pynanalogue.polars_bam_mods(str(bam_path))

    # Define the expected schema for polars_bam_mods output
    schema = {
        "read_id": pl.Utf8,
        "seq_len": pl.UInt64,
        "alignment_type": pl.Utf8,
        "align_start": pl.UInt64,
        "align_end": pl.UInt64,
        "contig": pl.Utf8,
        "contig_id": pl.Int32,
        "base": pl.Utf8,
        "is_strand_plus": pl.Boolean,
        "mod_code": pl.Utf8,
        "position": pl.UInt64,
        "ref_position": pl.Int64,
        "mod_quality": pl.UInt32,
    }

    # Read expected output as TSV with explicit schema
    expected_df = pl.read_csv(expected_path, separator="\t", schema=schema)

    # Compare DataFrames (works for both empty and non-empty DataFrames)
    assert result_df.equals(expected_df), (
        f"polars_bam_mods output doesn't match expected output for example_{example_number}"
    )


def test_example_1_bam_polars_bam_mods():
    """Test polars_bam_mods on example_1.bam matches expected TSV output"""
    _test_polars_bam_mods_against_expected(1)


def test_example_3_bam_polars_bam_mods():
    """Test polars_bam_mods on example_3.bam matches expected TSV output"""
    _test_polars_bam_mods_against_expected(3)


def test_example_7_bam_polars_bam_mods():
    """Test polars_bam_mods on example_7.bam matches expected TSV output"""
    _test_polars_bam_mods_against_expected(7)


def test_polars_bam_mods_basic(simple_bam):
    """Test polars_bam_mods returns valid Polars DataFrame"""
    df = pynanalogue.polars_bam_mods(str(simple_bam))

    assert isinstance(df, pl.DataFrame)
    assert len(df) > 0
    # Check all expected columns are present
    expected_columns = {
        "read_id",
        "seq_len",
        "alignment_type",
        "align_start",
        "align_end",
        "contig",
        "contig_id",
        "base",
        "is_strand_plus",
        "mod_code",
        "position",
        "ref_position",
        "mod_quality",
    }
    assert expected_columns.issubset(set(df.columns))


def test_polars_bam_mods_filters(simple_bam):
    """Test polars_bam_mods filtering works correctly,
    using one `InputBam` filter and one `InputMods` filter"""
    df_all = pynanalogue.polars_bam_mods(str(simple_bam), min_seq_len=0)
    df_filtered_1 = pynanalogue.polars_bam_mods(str(simple_bam), min_seq_len=6000)
    df_filtered_2 = pynanalogue.polars_bam_mods(str(simple_bam), mod_strand="bc_comp")

    # Our sequences are all length 5 kb, so a 6 kb filter should not let any through.
    assert len(df_filtered_1) == 0 and len(df_all) > 0

    # Mods are all on the basecalled strand, so we should get no information if we filter
    # to only allow the complement through
    assert len(df_filtered_2) == 0 and len(df_all) > 0


def test_polars_bam_mods_invalid_path():
    """Test proper error handling for missing files"""
    with pytest.raises(Exception):  # Update with specific exception
        pynanalogue.polars_bam_mods("/nonexistent/path.bam")
