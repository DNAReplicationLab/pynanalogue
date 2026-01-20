# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.6] - 2026-01-20

### Added
- New `peek()` function to extract BAM file metadata (contigs and modifications) without full processing
- New `seq_table()` function to extract read sequences and base qualities for a genomic region as a Polars DataFrame, with insertions shown as lowercase, deletions as periods, and modifications as 'Z'
- New `tag` parameter for `polars_bam_mods`, `read_info`, `window_reads`, and `seq_table` functions to filter by specific modification type (e.g., single-letter code "m" or ChEBI code "76792")
- New `win_op` parameter for `window_reads()` function supporting "density" (default) and "grad_density" modes to measure modification density or gradient within each window
- PyPI publishing via GitHub Actions using OIDC trusted publishing
- Comprehensive filtering parameter tests for `polars_bam_mods`, `read_info`, `window_reads`, and `seq_table` functions
- Gradient window tests with new example BAM files (`example_10`, `example_11`)
- New `two_mods_bam` test fixture for testing multiple modification types
- New pynanalogue-specific test data files for `seq_table` testing

### Changed
- Removed flaky 10% sampling test case in `seq_table` testing due to stochastic fluctuations with small sample sizes
- Updated nanalogue dependency from 0.1.4 to 0.1.6
- Consolidated `build-wheels.yml` into `publish_to_pypi.yml` workflow
- Re-enabled musl benchmark job in CI
- Renamed `test_filtering.py` to `test_polars_bam_mods_filtering.py` for clarity

### Fixed
- Fixed mismatch generation bug when simulating BAM files (via nanalogue 0.1.6)
- Fixed `full_region` misinterpretation when an open-ended interval or entire contig was provided in the `region` parameter (via nanalogue 0.1.6)

## [0.1.4] - 2026-01-11

### Added
- ARM64 (aarch64) support for Linux wheel builds across musllinux and manylinux platforms
- macOS x86\_64 wheel builds added

### Changed
- Updated nanalogue dependency from 0.1.2 to 0.1.4, allowing Mm/Ml tags in addition to MM/ML
- Removed Rust caching steps from CI/CD workflows to ensure clean builds
- Standard `cargo update` to update dependencies

## [0.1.0] - 2026-01-05

### Added
- Initial release of pynanalogue
- Python bindings for Nanalogue using PyO3
- Support for Python 3.9+
- Integration with polars for data manipulation
- SECURITY.md with security policy
- CONTRIBUTING.md with contribution guidelines
- CHANGELOG.md to track project changes
