# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
