# Tests for pynanalogue.window_reads() function
# Validates window-based modification analysis output

from io import StringIO
from pathlib import Path
import polars as pl
import pytest
import pynanalogue


def _test_window_reads_against_expected(example_number):
    """Helper function to test window_reads output against expected output for a given example

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
        / f"example_{example_number}_window_reads"
    )

    # Run window_reads with win=2, step=1
    result_df = pynanalogue.window_reads(str(bam_path), win=2, step=1)

    # Define the expected schema for window_reads output
    schema = {
        "contig": pl.Utf8,
        "ref_win_start": pl.Int64,
        "ref_win_end": pl.Int64,
        "read_id": pl.Utf8,
        "win_val": pl.Float32,
        "strand": pl.Utf8,
        "base": pl.Utf8,
        "mod_strand": pl.Utf8,
        "mod_type": pl.Utf8,
        "win_start": pl.UInt64,
        "win_end": pl.UInt64,
        "basecall_qual": pl.UInt32,
    }

    # Read expected output - first line has column names starting with '#'
    # Read the file and strip '#' from the first line before parsing
    with open(expected_path, "r") as f:
        lines = f.readlines()
        # Strip '#' from the first line (column names)
        lines[0] = lines[0].lstrip("#")

    # Create a temporary string buffer with the modified content
    modified_content = StringIO("".join(lines))

    # Read as TSV with explicit schema
    expected_df = pl.read_csv(modified_content, separator="\t", schema=schema)

    # Compare DataFrames (sort both to ensure consistent ordering)
    result_sorted = result_df.sort(by=result_df.columns)
    expected_sorted = expected_df.sort(by=expected_df.columns)

    assert result_sorted.equals(expected_sorted), (
        "DataFrame contents don't match expected output"
    )


def test_example_1_bam_window_reads():
    """Test window_reads on example_1.bam matches expected output"""
    _test_window_reads_against_expected(1)


def test_example_3_bam_window_reads():
    """Test window_reads on example_3.bam matches expected output"""
    _test_window_reads_against_expected(3)


def test_example_7_bam_window_reads():
    """Test window_reads on example_7.bam matches expected output"""
    _test_window_reads_against_expected(7)


def _test_window_reads_gradient_against_expected(example_number, win, step):
    """Helper function to test window_reads gradient output against expected output

    Args:
        example_number: The example number (e.g., 10, 11)
        win: Window size parameter
        step: Step size parameter
    """
    # Paths to test data
    test_dir = Path(__file__).parent
    bam_path = test_dir / "data" / "examples" / f"example_{example_number}.bam"
    expected_path = (
        test_dir
        / "data"
        / "expected_outputs"
        / f"example_{example_number}_win_grad_win_{win}_step_{step}"
    )

    # Run window_reads with win_op="grad_density"
    result_df = pynanalogue.window_reads(
        str(bam_path), win=win, step=step, win_op="grad_density"
    )

    # Define the expected schema for window_reads output
    schema = {
        "contig": pl.Utf8,
        "ref_win_start": pl.Int64,
        "ref_win_end": pl.Int64,
        "read_id": pl.Utf8,
        "win_val": pl.Float32,
        "strand": pl.Utf8,
        "base": pl.Utf8,
        "mod_strand": pl.Utf8,
        "mod_type": pl.Utf8,
        "win_start": pl.UInt64,
        "win_end": pl.UInt64,
        "basecall_qual": pl.UInt32,
    }

    # Read expected output - first line has column names starting with '#'
    with open(expected_path, "r") as f:
        lines = f.readlines()
        # Strip '#' from the first line (column names)
        lines[0] = lines[0].lstrip("#")

    # Create a temporary string buffer with the modified content
    modified_content = StringIO("".join(lines))

    # Read as TSV with explicit schema
    expected_df = pl.read_csv(modified_content, separator="\t", schema=schema)

    # Compare DataFrames (sort both to ensure consistent ordering)
    result_sorted = result_df.sort(by=result_df.columns)
    expected_sorted = expected_df.sort(by=expected_df.columns)

    assert result_sorted.equals(expected_sorted), (
        "DataFrame contents don't match expected output"
    )


def test_example_10_bam_window_reads_gradient_win10_step1():
    """Test window_reads gradient on example_10.bam with win=10, step=1"""
    _test_window_reads_gradient_against_expected(10, win=10, step=1)


def test_example_10_bam_window_reads_gradient_win20_step2():
    """Test window_reads gradient on example_10.bam with win=20, step=2"""
    _test_window_reads_gradient_against_expected(10, win=20, step=2)


def test_example_11_bam_window_reads_gradient_win10_step1():
    """Test window_reads gradient on example_11.bam with win=10, step=1"""
    _test_window_reads_gradient_against_expected(11, win=10, step=1)


def test_example_11_bam_window_reads_gradient_win20_step2():
    """Test window_reads gradient on example_11.bam with win=20, step=2"""
    _test_window_reads_gradient_against_expected(11, win=20, step=2)


def test_window_reads_invalid_win_op_raises_error(simple_bam):
    """Test that an invalid win_op value raises an error"""
    with pytest.raises(
        Exception, match="win_op must be set to density or grad_density"
    ):
        pynanalogue.window_reads(
            str(simple_bam), win=5, step=2, win_op="invalid_option"
        )


def test_window_reads_basic(simple_bam):
    """Test window_reads returns valid Polars DataFrame"""
    df = pynanalogue.window_reads(str(simple_bam), win=5, step=2)

    assert isinstance(df, pl.DataFrame)
    assert len(df) >= 0  # Could be 0 if no mods, but should be valid

    # Check for expected columns based on the example output in lib.rs
    expected_cols = [
        "contig",
        "ref_win_start",
        "ref_win_end",
        "read_id",
        "win_val",
        "strand",
        "base",
        "mod_strand",
        "mod_type",
    ]

    for col in expected_cols:
        assert col in df.columns, f"Expected column '{col}' not found"
