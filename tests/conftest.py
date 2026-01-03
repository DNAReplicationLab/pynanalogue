# Pytest configuration and shared fixtures for pynanalogue tests
# Provides reusable test fixtures including BAM file generation

import pytest
import json
from pathlib import Path
from uuid import uuid4
import pynanalogue


def _generate_bam_from_config(tmp_path, config_name, prefix="test"):
    """Helper function to generate BAM files from simulation configs

    Args:
        tmp_path: pytest tmp_path fixture
        config_name: Name of the config file (without .json extension)
        prefix: Prefix for generated BAM/FASTA filenames

    Returns:
        Path to the generated BAM file
    """
    # Load config from file
    config_file = (
        Path(__file__).parent / "data" / "simulation_configs" / f"{config_name}.json"
    )

    # Check if config file exists
    if not config_file.exists():
        raise RuntimeError(
            f"Simulation config file not found: {config_file}\n"
            f"Expected path: {config_file.absolute()}"
        )

    config = json.loads(config_file.read_text())

    # Use unique filenames for concurrent test execution
    unique_id = uuid4().hex[:8]
    bam_path = tmp_path / f"{prefix}_{unique_id}.bam"
    fasta_path = tmp_path / f"{prefix}_{unique_id}.fasta"

    pynanalogue.simulate_mod_bam(
        json_config=json.dumps(config),
        bam_path=str(bam_path),
        fasta_path=str(fasta_path),
    )

    return bam_path


@pytest.fixture
def simple_bam(tmp_path):
    """Generate a simple test BAM using `simulate_mod_bam`"""
    return _generate_bam_from_config(tmp_path, "simple_bam", "test")


@pytest.fixture
def benchmark_bam(tmp_path):
    """Generate a larger benchmark BAM for performance testing"""
    return _generate_bam_from_config(tmp_path, "benchmark_bam", "benchmark")
