# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.6] - 2026-01-19

### Added
- New `peek()` function to extract BAM file metadata (contigs and modifications) without full processing
- New `tag` parameter for `polars_bam_mods`, `read_info`, and `window_reads` functions to filter by specific modification type (e.g., single-letter code "m" or ChEBI code "76792")
- Comprehensive filtering parameter tests for `polars_bam_mods`, `read_info`, and `window_reads` functions
- New `two_mods_bam` test fixture for testing multiple modification types

### Changed
- Updated nanalogue dependency from 0.1.4 to 0.1.6
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
