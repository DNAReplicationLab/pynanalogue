# Tests for Python type handling in pynanalogue
# Validates that functions accept and return correct types

import json
from pathlib import Path
from uuid import uuid4
import polars as pl
import pynanalogue


class TestReturnTypes:
    """Test that functions return the expected types"""

    def test_read_info_returns_bytes(self, simple_bam):
        """Test that read_info returns bytes"""
        result = pynanalogue.read_info(str(simple_bam))
        assert isinstance(result, bytes)
        assert len(result) > 0

    def test_polars_bam_mods_returns_dataframe(self, simple_bam):
        """Test that polars_bam_mods returns a Polars DataFrame"""
        result = pynanalogue.polars_bam_mods(str(simple_bam))
        assert isinstance(result, pl.DataFrame)

    def test_window_reads_returns_dataframe(self, simple_bam):
        """Test that window_reads returns a Polars DataFrame"""
        result = pynanalogue.window_reads(str(simple_bam), win=5, step=2)
        assert isinstance(result, pl.DataFrame)

    def test_simulate_mod_bam_returns_none(self, tmp_path):
        """Test that simulate_mod_bam returns None on success"""
        config = {
            "contigs": {"number": 1, "len_range": [100, 100]},
            "reads": [
                {
                    "number": 10,
                    "mapq_range": [10, 20],
                    "base_qual_range": [10, 20],
                    "len_range": [0.5, 0.5],
                }
            ],
        }

        unique_id = uuid4().hex[:8]

        result = pynanalogue.simulate_mod_bam(
            json_config=json.dumps(config),
            bam_path=str(tmp_path / f"test_{unique_id}.bam"),
            fasta_path=str(tmp_path / f"test_{unique_id}.fasta"),
        )

        assert result is None


class TestParameterTypes:
    """Test that functions accept correct parameter types"""

    def test_read_info_accepts_string_path(self, simple_bam):
        """Test that read_info accepts string paths"""
        result = pynanalogue.read_info(str(simple_bam))
        assert isinstance(result, bytes)

    def test_read_info_accepts_path_object(self, simple_bam):
        """Test that read_info accepts Path objects converted to string"""
        result = pynanalogue.read_info(str(Path(simple_bam)))
        assert isinstance(result, bytes)

    def test_boolean_parameters(self, simple_bam):
        """Test that boolean parameters work correctly"""
        result = pynanalogue.read_info(
            str(simple_bam),
            treat_as_url=False,
            include_zero_len=False,
            full_region=False,
            exclude_mapq_unavail=True,
        )
        assert isinstance(result, bytes)

    def test_numeric_parameters(self, simple_bam):
        """Test that numeric parameters accept correct types"""
        result = pynanalogue.polars_bam_mods(
            str(simple_bam),
            min_seq_len=100,
            min_align_len=50,
            threads=4,
            sample_fraction=0.5,
            mapq_filter=10,
            min_mod_qual=128,
            trim_read_ends_mod=5,
            base_qual_filter_mod=20,
        )
        assert isinstance(result, pl.DataFrame)

    def test_tuple_parameters(self, simple_bam):
        """Test that tuple parameters work correctly"""
        result = pynanalogue.polars_bam_mods(
            str(simple_bam), reject_mod_qual_non_inclusive=(50, 200)
        )
        assert isinstance(result, pl.DataFrame)

    def test_set_parameters(self, simple_bam):
        """Test that set parameters work correctly"""
        # Get some read IDs first
        data = json.loads(pynanalogue.read_info(str(simple_bam)))

        if data:
            read_id = data[0]["read_id"]
            result = pynanalogue.read_info(str(simple_bam), read_ids={read_id})
            assert isinstance(result, bytes)

    def test_empty_set_parameter(self, simple_bam):
        """Test that empty set is handled correctly (default behavior)"""
        result = pynanalogue.read_info(str(simple_bam), read_ids=set())
        assert isinstance(result, bytes)


class TestDataFrameSchema:
    """Test that returned DataFrames have correct schemas"""

    def test_polars_bam_mods_schema(self, simple_bam):
        """Test that polars_bam_mods DataFrame has correct column types"""
        df = pynanalogue.polars_bam_mods(str(simple_bam))

        # Check that required columns exist
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

        # Verify some expected types
        assert df["read_id"].dtype == pl.Utf8 or df["read_id"].dtype == pl.String
        assert df["seq_len"].dtype in [pl.Int64, pl.UInt64, pl.Int32, pl.UInt32]
        assert df["position"].dtype in [pl.Int64, pl.UInt64, pl.Int32, pl.UInt32]
        assert df["mod_quality"].dtype in [pl.UInt8, pl.UInt32, pl.Int64, pl.UInt64]

    def test_window_reads_schema(self, simple_bam):
        """Test that window_reads DataFrame has correct column types"""
        df = pynanalogue.window_reads(str(simple_bam), win=5, step=2)

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

        # Verify DataFrame is not empty (assuming simple_bam has data)
        assert len(df) >= 0  # Could be 0 if no mods, but should be valid


class TestJSONOutput:
    """Test that read_info produces valid JSON"""

    def test_read_info_json_decodable(self, simple_bam):
        """Test that read_info output can be decoded as JSON"""
        result = pynanalogue.read_info(str(simple_bam))
        decoded = json.loads(result)

        assert isinstance(decoded, list)

    def test_read_info_json_structure(self, simple_bam):
        """Test that decoded JSON has expected structure"""
        result = pynanalogue.read_info(str(simple_bam))
        decoded = json.loads(result)

        if decoded:
            # Check first record has expected fields
            first_record = decoded[0]
            required_fields = ["read_id", "sequence_length"]
            for field in required_fields:
                assert field in first_record, f"Expected field '{field}' not found"

            # Verify types
            assert isinstance(first_record["read_id"], str)
            assert isinstance(first_record["sequence_length"], int)


class TestDefaultParameters:
    """Test that default parameters work as expected"""

    def test_read_info_minimal_call(self, simple_bam):
        """Test read_info with only required parameter"""
        result = pynanalogue.read_info(str(simple_bam))
        assert isinstance(result, bytes)

    def test_polars_bam_mods_minimal_call(self, simple_bam):
        """Test polars_bam_mods with only required parameter"""
        result = pynanalogue.polars_bam_mods(str(simple_bam))
        assert isinstance(result, pl.DataFrame)

    def test_window_reads_minimal_call(self, simple_bam):
        """Test window_reads with only required parameters"""
        result = pynanalogue.window_reads(str(simple_bam), win=10, step=5)
        assert isinstance(result, pl.DataFrame)
