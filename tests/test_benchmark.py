# Performance benchmarks for pynanalogue functions using pytest-benchmark
# Tests use benchmark_bam.json config for reproducible results across hardware

import pynanalogue


class TestBenchmarks:
    """Performance benchmarks for pynanalogue functions"""

    def test_benchmark_polars_bam_mods(self, benchmark, benchmark_bam):
        """Benchmark polars_bam_mods on 5k reads, 1 Mbp contigs"""
        # The benchmark fixture automatically times the function call
        result = benchmark(pynanalogue.polars_bam_mods, str(benchmark_bam))
        assert len(result) > 0

    def test_benchmark_read_info(self, benchmark, benchmark_bam):
        """Benchmark read_info on 5k reads"""
        result = benchmark(pynanalogue.read_info, str(benchmark_bam))
        assert len(result) > 0

    def test_benchmark_window_reads(self, benchmark, benchmark_bam):
        """Benchmark window_reads on 5k reads with 5bp window"""
        result = benchmark(pynanalogue.window_reads, str(benchmark_bam), win=5, step=2)
        assert len(result) > 0
