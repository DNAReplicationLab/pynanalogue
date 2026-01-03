# Tests for error handling in pynanalogue Python bindings
# Focuses on PyO3 interface errors, not business logic errors

import pytest
import pynanalogue

# This test file validates Python-Rust binding error handling -
# focusing only on errors that occur at the PyO3 interface layer,
# not business logic errors from the underlying `nanalogue` library.
# We expect `nanalogue` is robust enough to deal with errors such
# as malformed regions or numbers in the wrong order etc.


class TestPythonTypeErrors:
    """Test that passing wrong Python types raises appropriate errors"""

    def test_bam_path_not_string(self):
        """Test that passing non-string bam_path raises TypeError"""
        with pytest.raises(TypeError):
            pynanalogue.read_info(12345)  # int instead of str

        with pytest.raises(TypeError):
            pynanalogue.read_info(None)  # None instead of str

    def test_boolean_params_not_bool(self):
        """Test that boolean parameters reject non-boolean values"""
        with pytest.raises(TypeError):
            pynanalogue.read_info(
                "/some/path.bam",
                treat_as_url="true",  # string instead of bool
            )

    def test_numeric_params_wrong_type(self):
        """Test that numeric parameters reject wrong types"""
        with pytest.raises(TypeError):
            pynanalogue.read_info(
                "/some/path.bam",
                min_seq_len="100",  # string instead of int
            )

        with pytest.raises(TypeError):
            pynanalogue.polars_bam_mods(
                "/some/path.bam",
                sample_fraction="0.5",  # string instead of float
            )

    def test_tuple_param_wrong_type(self):
        """Test that tuple parameters reject wrong types"""
        with pytest.raises(TypeError):
            pynanalogue.polars_bam_mods(
                "/some/path.bam",
                reject_mod_qual_non_inclusive=[50, 200],  # list instead of tuple
            )

    def test_set_param_wrong_type(self):
        """Test that set parameters reject wrong types"""
        with pytest.raises(TypeError):
            pynanalogue.read_info(
                "/some/path.bam",
                read_ids=["id1", "id2"],  # list instead of set
            )


class TestWindowReadsRequiredParams:
    """Test that window_reads enforces required parameters"""

    def test_missing_win_parameter(self):
        """Test that win parameter is required"""
        with pytest.raises(TypeError):
            pynanalogue.window_reads("/some/path.bam", step=5)

    def test_missing_step_parameter(self):
        """Test that step parameter is required"""
        with pytest.raises(TypeError):
            pynanalogue.window_reads("/some/path.bam", win=10)


class TestSimulateModBamRequiredParams:
    """Test that simulate_mod_bam enforces required parameters"""

    def test_missing_json_config(self):
        """Test that json_config is required"""
        with pytest.raises(TypeError):
            pynanalogue.simulate_mod_bam(
                bam_path="/tmp/test.bam", fasta_path="/tmp/test.fasta"
            )

    def test_missing_bam_path(self):
        """Test that bam_path is required"""
        with pytest.raises(TypeError):
            pynanalogue.simulate_mod_bam(json_config="{}", fasta_path="/tmp/test.fasta")

    def test_missing_fasta_path(self):
        """Test that fasta_path is required"""
        with pytest.raises(TypeError):
            pynanalogue.simulate_mod_bam(json_config="{}", bam_path="/tmp/test.bam")

    def test_json_config_not_string(self):
        """Test that json_config must be a string"""
        with pytest.raises(TypeError):
            pynanalogue.simulate_mod_bam(
                json_config={"key": "value"},  # dict instead of string
                bam_path="/tmp/test.bam",
                fasta_path="/tmp/test.fasta",
            )


class TestZeroLengthReadGuard:
    """Test that include_zero_len parameter is properly guarded"""

    def test_include_zero_len_raises_error(self):
        """Test that include_zero_len=True raises ValueError with clear message"""
        with pytest.raises(
            ValueError,
            match="include_zero_len=True is not yet supported due to potential crashes",
        ):
            pynanalogue.read_info("/some/path.bam", include_zero_len=True)

        with pytest.raises(
            ValueError,
            match="include_zero_len=True is not yet supported due to potential crashes",
        ):
            pynanalogue.polars_bam_mods("/some/path.bam", include_zero_len=True)

        with pytest.raises(
            ValueError,
            match="include_zero_len=True is not yet supported due to potential crashes",
        ):
            pynanalogue.window_reads(
                "/some/path.bam", win=5, step=2, include_zero_len=True
            )
