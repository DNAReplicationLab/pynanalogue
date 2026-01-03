# Contributing to pynanalogue

Thank you for your interest in contributing to pynanalogue! This document provides guidelines for contributing to the project.

## Code of Conduct

Please be respectful and constructive in all interactions with other contributors and maintainers.

## Development Setup

### Prerequisites

- Rust 1.70 or higher
- Python 3.10 or higher
- maturin for building Python wheels

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
   pip install -e .[dev]  # Once optional dependencies are added
   ```

4. **Build the project in development mode**:
   ```bash
   maturin develop
   ```

5. **Verify the installation**:
   ```bash
   python -c "import pynanalogue; print('Success!')"
   ```

## Making Changes

### Code Style

#### Rust Code
- Follow standard Rust formatting: run `cargo fmt` before committing
- Pass all clippy lints: `cargo clippy --all-targets --all-features -- -D warnings`
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
cargo test
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

### Building the Project

**Development build** (faster, with debug symbols):
```bash
maturin develop
```

**Release build** (optimized):
```bash
maturin build --release
```

**Check for issues without building**:
```bash
cargo check
cargo clippy
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
   cargo clippy --all-targets --all-features -- -D warnings
   cargo test

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

When adding development dependencies, use optional dependency groups:

```toml
[project.optional-dependencies]
dev = [
    "pytest>=7.0",
    "pytest-cov",
    "pytest-benchmark",
    "ruff>=0.1.0",
]
```

Install with: `pip install -e .[dev]`

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