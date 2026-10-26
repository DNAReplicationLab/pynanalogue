# Contributing to pynanalogue

Thank you for your interest in contributing to pynanalogue! This document provides guidelines for contributing to the project.

## Code of Conduct

Please be respectful and constructive in all interactions with other contributors and maintainers.

## Development Setup

### Prerequisites

- Current stable Rust (edition 2024; a minimum supported Rust version for the full dependency stack has not been established)
- Python 3.11 or higher
- maturin for building Python wheels
- Native build tools: a C compiler, clang/libclang for bindgen, pkg-config, and Perl
- Development headers/libraries for compression (zlib, bzip2, xz), OpenSSL, and libcurl; package names vary by operating system

### Setting Up Your Development Environment

1. **Clone the repository**:
   ```bash
   git clone https://github.com/DNAReplicationLab/pynanalogue.git
   cd pynanalogue
   ```

2. **Create and activate a virtual environment**:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```

3. **Install development dependencies**:
   ```bash
   pip install maturin
   pip install -e '.[dev]' --config-settings=build-args=--locked
   ```

4. **Build the project in development mode**:
   ```bash
   maturin develop --locked
   ```

5. **Verify the installation**:
   ```bash
   python -c "import pynanalogue; print('Success!')"
   ```

## Making Changes

### Code Style

#### Rust Code
- Follow standard Rust formatting: run `cargo fmt` before committing
- Pass all clippy lints: `cargo clippy --locked --all-targets --all-features -- -D warnings`
- This project uses 80+ strict clippy lints - they're intentional and helpful for code quality
- All code files should start with a brief 2-line comment explaining what the file does
- Add proper documentation to public APIs

#### Python Code
- Follow PEP 8 guidelines
- Use `ruff` for linting and formatting
- Run `ruff check tests/` to check for linting issues
- Run `ruff format tests/` to format code
- Maintain consistency with existing code style in each file

### Testing

**Tests are required for all new functionality and bug fixes.**

This project follows Test-Driven Development (TDD) principles:

1. **Write tests before implementing functionality**
2. **Write minimal code to make the test pass**
3. **Refactor while keeping tests green**
4. **Ensure all tests pass before submitting a PR**

#### Running Tests

**Rust tests**:
```bash
cargo test --locked
```

**Python tests**:
```bash
pytest tests/
pytest --cov=pynanalogue tests/  # With coverage
```

**Python linting**:
```bash
ruff check tests/              # Check for linting issues
ruff check tests/ --fix        # Auto-fix linting issues
ruff format tests/             # Format code
ruff format --check tests/     # Check formatting without changing files
```

#### Test Structure
- Add test data to the `tests/data/` directory as needed
- Unit tests should cover individual functions and edge cases
- Integration tests should verify end-to-end functionality

#### Smoke-testing release wheels

Run **Smoke-test built wheels** manually in GitHub Actions and enter the run ID
of a successful **Build Wheels and Publish to PyPI** run (including build-only
runs). The artifacts must still be available.

This installs the existing wheels in fresh environments on native x86_64 and
ARM64 Linux and macOS runners, using matching manylinux/musllinux containers on
Linux. It checks the installed version against the build revision and runs that
revision's functional tests, including all six APIs and simulation/read-back,
on Python 3.11 (the stable ABI baseline). Benchmarks are excluded. Runtime
dependencies are resolved from PyPI as binary wheels; this is not a minimum
dependency-version test.

The workflow does not rebuild, publish, or change tags. It is a separate manual
check, not yet a publishing gate; it does not test the source distribution or
every supported Python version.

### Building the Project

**Development build** (unoptimized; debug symbols are disabled by the repository's dev profile):
```bash
maturin develop --locked
```

**Release build** (optimized):
```bash
maturin build --release --locked
```

**Check for issues without building**:
```bash
cargo check --locked
cargo clippy --locked
```

### Commit Messages

- Use clear, descriptive commit messages
- Start with a verb such as "adds", "fixes", "updates" etc.
- Keep the first line under 50 characters
- Add a blank line, then detailed explanation if needed
- Reference issue numbers when applicable (e.g., "Fixes #123")

Example:
```
Add support for parsing ModBAM files

Implements new parsing logic to handle modification tags
in BAM files. Includes tests for various mod types.

Fixes #42
```

### Pull Request Process

1. **Fork the repository** and create a new branch for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** following the code style guidelines

3. **Add or update tests** as needed - all tests must pass

4. **Run code quality checks**:
   ```bash
   # Rust checks
   cargo fmt
   cargo clippy --locked --all-targets --all-features -- -D warnings
   cargo test --locked

   # Python checks
   ruff check tests/
   ruff format tests/
   pytest tests/
   ```

5. **Update documentation** if you're changing functionality:
   - Update README.md if needed
   - Update CHANGELOG.md following Keep a Changelog format
   - Add docstrings/doc comments to new code

6. **Commit your changes** with clear commit messages

7. **Push to your fork** and submit a pull request

8. **Address review feedback** - maintainers may request changes

### Pre-commit Hooks

If pre-commit hooks are configured:
- **Ensure all hooks pass before committing**
- **Never use `--no-verify`** to bypass hooks
- If hooks fail, fix the issues they identify
- Bypassing quality checks is not acceptable

## Reporting Issues

When reporting issues, please include:

- **Clear description** of the problem
- **Steps to reproduce** the issue
- **Expected behavior** vs **actual behavior**
- **Environment details**:
  - Operating System
  - Python version (`python --version`)
  - Rust version (`rustc --version`)
  - Package version
- **Error messages or logs** (full stack traces help!)
- **Test data or examples** that reproduce the issue (if applicable)

## Feature Requests

We welcome feature requests! Please open an issue describing:

- **The problem** you're trying to solve
- **Your proposed solution** and how it would work
- **Alternatives you've considered**
- **Use cases** - how would this feature be used?
- **Willingness to implement** - can you contribute the implementation?

## Documentation

Good documentation is essential:

- **Update README.md** when adding or changing functionality
- **Add docstrings/doc comments** to new functions and classes
- **Update CHANGELOG.md** following Keep a Changelog format
- **Add examples** to demonstrate usage when appropriate
- **Keep documentation current** - outdated docs are worse than no docs

## Development Dependencies

The existing `dev` optional dependency group in `pyproject.toml` includes pytest,
pytest-cov, pytest-benchmark, and ruff. Add development dependencies there.

Install with: `pip install -e '.[dev]' --config-settings=build-args=--locked`

## Design Philosophy

This project prioritizes **Python ergonomics over Rust elegance**:

- Functions may have many parameters for explicit, clear Python APIs
- Some Rust code patterns are intentionally avoided for better Python UX
- Error messages should be helpful to Python users, not just Rust developers
- See module-level documentation in `src/lib.rs` for detailed rationale

## Versioning

This project follows [Semantic Versioning](https://semver.org/):

- **Pre-1.0 (current)**: API may change without notice in minor versions
- **Post-1.0**: Breaking changes only in major versions

## Questions?

If you have questions about contributing:
- Open an issue with the "question" label
- Check existing issues and pull requests
- Review the README and existing documentation

## License

By contributing to pynanalogue, you agree that your contributions will be licensed under the MIT License, the same license as the project.

---

Thank you for contributing to pynanalogue! Your efforts help make bioinformatics tools better for everyone.
