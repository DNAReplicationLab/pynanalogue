# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.1] - 2026-10-26

**Upgrading from 0.1.7:** Python 3.11 or newer is now required. Several DataFrame
columns use narrower integer types, and input validation is stricter; see the
compatibility details under Changed below. Python function signatures are unchanged.

### Added
- Mapping quality (`mapq`) in `read_info()` JSON output and `polars_bam_mods()` DataFrames, with updated documentation and examples
- `simulate_mod_bam()` accepts a JSON `seed` for reproducible simulation and `mm_suffix`/`drop` modification options; a `.cram` output filename selects CRAM 3.1 output with a generated reference and index
- Default-on build-only option for manual release workflow runs (CI-only; no package-code changes), retaining wheels and source archives as artifacts without publishing to PyPI
- Manual wheel smoke-test workflow (CI-only; no package-code changes) that tests existing build artifacts on native Linux/macOS architectures without rebuilding or publishing
- Regression tests for numeric input limits, MAPQ filtering, exact output schemas (including empty results), and propagation of analysis errors
- Dependency-security CI on pull requests, main-branch pushes, weekly schedules, and manual runs using pinned cargo-audit and cargo-deny versions (advisories, bans, and sources; no license policy)
- Explicit advisory exceptions for unmaintained custom_derive in both security tools, and for inactive bincode and quick-xml dependencies in cargo-audit, with documented reasons and conditions for revisiting them
- Scoped exceptions for PyO3 advisories RUSTSEC-2026-0176 and RUSTSEC-2026-0177 in both security tools, backed by a documented reachability investigation in SECURITY.md

### Changed
- Raised the minimum Python version to 3.11, updated the stable ABI baseline to Python 3.11, and moved CI and wheel-build interpreter coverage to Python 3.11–3.15
- Changed Polars output types from `UInt64` to `UInt32` for `window_reads()` columns `win_start` and `win_end`, and `polars_bam_mods()` columns `seq_len`, `align_start`, `align_end`, and `position`
- Limited `min_seq_len`, `min_align_len`, `win`, and `step` to at most 4294967295 (u32::MAX); larger values now raise `ValueError`
- Inherited stricter BAM/modification and identifier validation, including maximum read length 16,777,215 bases, maximum supported contig length 4,294,967,294, and at most 20 modification groups per record
- Region strings reject explicit numeric ends above 4,294,967,294; use an omitted end such as `chr1:0-` for an open-ended region. `seq_table()` also rejects oversized header contig lengths when resolving a region
- `peek()` now raises an error for input containing no records; `read_info()` still returns an empty JSON array, while the table/window APIs return typed empty DataFrames for valid no-result cases
- Simulator output filenames must end case-insensitively in `.bam` or `.cram`; simulated unmapped reads now have MAPQ 0 rather than 255
- URL readers now consistently restrict schemes to HTTP, HTTPS and FTP
- Optimized modification parsing and overlapping-window calculations; window `basecall_qual` can differ by one Phred unit at rounding boundaries because floating-point summation order changed
- Updated nanalogue from 0.1.7 through 0.1.11 and 0.2.0 to an exact 0.2.1 dependency with its `polars` feature enabled, refreshed Cargo.lock, and reassessed the existing advisory exceptions
- Updated Rust dependencies: polars to 0.53, rust-htslib to 1.0.0, pyo3 to 0.27, and pyo3-polars to 0.26.0
- Aligned MAPQ filtering with nanalogue 0.2.0: MAPQ 0 is subject to the numeric threshold, while MAPQ 255 denotes unavailable quality and is excluded only when `exclude_mapq_unavail` is enabled
- Disabled debug symbols for development builds and their dependencies, and made the test profile inherit the development profile
- Pinned GitHub Actions to full commit hashes, restricted workflow token permissions, and disabled checkout credential persistence
- Enforced Cargo.lock for normal Rust checks and maturin compilation in tests, benchmarks, and wheel builds; retained the separate fresh-dependency Rust workflow
- Added a seven-day cooldown for Dependabot version updates across Cargo, Python, and GitHub Actions
- Hardened CI shell interpolation, pinned Linux container images by digest, documented publishing permissions, replaced redundant Rust setup actions, and limited concurrent runs while preserving active publications

### Fixed
- Build Linux release wheels on Ubuntu using maturin's cargo-zigbuild integration, with matching clang/libclang for HTSlib bindings and the x86_64 libdeflate compatibility flags; retain the existing manylinux and musllinux targets
- Reject malformed MM/ML/MN data instead of panicking or silently accepting inconsistent records; handle canonical `U` as BAM-encoded `T`, support padding (`P`) CIGAR operations, and validate explicitly empty MM tags
- Fixed `mod_region` filtering when the requested region contains or exactly matches the whole read alignment: terminal soft-clipped modification calls are excluded, while insertions between retained calls are still kept
- Reject invalid simulation schedules and aliased output paths before writing output
- Fixed HTTPS certificate discovery on macOS. HTTPS access still requires an up-to-date system or custom CA bundle; upgrading pynanalogue does not refresh that trust store
- Updated compatible transitive dependencies to resolve security advisories and replace yanked chacha20 0.10.1; remaining PyO3 and quick-xml advisories have documented exceptions rather than dependency fixes
- Fixed builds on Rust 1.99 by updating ethnum to 1.5.3, and resolved Ruff findings in test imports, set annotations, and the missing-file exception assertion
- Return typed, zero-row DataFrames from `window_reads()` and `seq_table()` when filtering leaves no records, or windowing yields no rows; unrelated processing errors still propagate
- Fixed musl benchmark CI builds by installing Alpine-compatible Clang 20 and its development libraries, LLVM, and Perl, and configuring Clang/libclang paths

## [0.1.7] - 2026-02-09

### Changed
- Removed redundant whitespace trimming for read_filter in Python bindings (handled by nanalogue library)
- Updated README code examples to show complete, reproducible output matching actual program behaviour
- Bumps nanalogue version to v0.1.7, this fixes a peek bug when zero len reads are encountered.

### Added
- Version check job in PyPI publish workflow (internal GitHub Actions change) to verify tag, Cargo.toml, Cargo.lock, and pyproject.toml versions match before building wheels
- Automated README example testing (`tests/test_readme_examples.py`) that extracts code blocks and expected output from HTML comment markers in README.md and verifies they match

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
