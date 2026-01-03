# Tests for pynanalogue.read_info() function
# Validates JSON output, filtering, and integration with example data

import json
from pathlib import Path
import pynanalogue


def _test_read_info_against_expected(example_number):
    """Helper function to test read_info output against expected JSON output for a given example

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
        / f"example_{example_number}_read_info.json"
    )

    # Run read_info
    result_bytes = pynanalogue.read_info(str(bam_path))
    result_json = json.loads(result_bytes)

    # Read expected output
    with open(expected_path, "r") as f:
        expected_json = json.load(f)

    # Compare JSON outputs
    assert result_json == expected_json, (
        f"read_info output doesn't match expected output for example_{example_number}"
    )


def test_example_1_bam_read_info():
    """Test read_info on example_1.bam matches expected JSON output"""
    _test_read_info_against_expected(1)


def test_example_3_bam_read_info():
    """Test read_info on example_3.bam matches expected JSON output"""
    _test_read_info_against_expected(3)


def test_example_7_bam_read_info():
    """Test read_info on example_7.bam matches expected JSON output"""
    _test_read_info_against_expected(7)


def test_read_info_basic(simple_bam):
    """Test read_info returns valid JSON bytes"""
    result = pynanalogue.read_info(str(simple_bam))

    assert isinstance(result, bytes)
    assert len(result) > 0

    # Verify it's valid JSON
    decoded = json.loads(result)
    assert isinstance(decoded, list)


def test_read_info_with_filters(simple_bam):
    """Test read_info with various filters"""
    # Test with min_seq_len filter
    result_all = pynanalogue.read_info(str(simple_bam))
    result_filtered = pynanalogue.read_info(str(simple_bam), min_seq_len=6000)
    result_filtered_2 = pynanalogue.read_info(str(simple_bam), mod_strand="bc_comp")

    decoded_all = json.loads(result_all)
    decoded_filtered = json.loads(result_filtered)
    decoded_filtered_2 = json.loads(result_filtered_2)

    # With 6kb filter, should have zero reads
    assert len(decoded_filtered) == 0 and len(decoded_all) > 0

    # Check mod_count fields in decoded_filtered_2 (with mod_strand="bc_comp")
    # All mod_count fields should be "NA"
    assert all(read["mod_count"] == "NA" for read in decoded_filtered_2)

    # Check mod_count fields in decoded_all
    # All mod_count fields should start with "T+T:" followed by a number
    assert all(read["mod_count"].startswith("T+T:") for read in decoded_all)
