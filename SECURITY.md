# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 0.x     | :x: Alpha          |

**Note**: This project is currently in alpha development. 
No versions are officially supported at this time.
Use at your own risk.

## Reporting a Vulnerability

If you discover a security vulnerability in `pynanalogue`, please report it responsibly:

1. **Do NOT open a public issue** for security vulnerabilities
2. Email the maintainer directly at: mail AT unintegrable dot com
3. Include:
   - Description of the vulnerability
   - Steps to reproduce
   - Potential impact
   - Any suggested fixes (optional)

We aim to respond to security reports as soon as possible.

## Security Considerations

### Data Sensitivity

This tool processes BAM files which may contain:
- Sensitive genetic/genomic data
- Personal health information
- Research data subject to confidentiality agreements

**Users are responsible for**:
- Ensuring proper access controls on input/output files
- Complying with data protection regulations (GDPR, HIPAA, etc.)
- Maintaining appropriate file permissions in multi-user environments

### URL Loading Feature

When using the `treat_as_url=True` option:
- Only use trusted URLs and secure (HTTPS) connections when possible
- Be aware that BAM files and indices will be downloaded/accessed over the network
- Ensure your network environment is secure
- Consider the privacy implications of accessing data from remote servers

### Dependencies

This package relies on:
- `nanalogue` (Rust library) - for core BAM processing
- `rust-htslib` - for BAM file I/O
- `polars` - for data frame operations

For a complete list of dependencies, see `Cargo.toml`.
Security of these dependencies is maintained upstream.
Users should keep dependencies updated.

### PyO3 advisory exceptions

On 2026-10-22, a [reachability investigation](https://ampcode.com/threads/T-01a1213c-352d-752f-b3dd-0063f9a9957a)
found no path through our six public Python functions to the operations affected by
[RUSTSEC-2026-0176](https://rustsec.org/advisories/RUSTSEC-2026-0176.html)
and [RUSTSEC-2026-0177](https://rustsec.org/advisories/RUSTSEC-2026-0177.html).
These two advisory IDs are therefore explicitly ignored in `deny.toml` and
`.cargo/audit.toml`. This is a non-reachability exception, not a dependency fix.

The reviewed stack was nanalogue 0.2.0, Rust Polars 0.53.0,
pyo3-polars 0.26.0, and PyO3 0.27.2 with `abi3-py311` and
`extension-module`. The investigation checked the relevant published crate
archives against Cargo.lock checksums, verified their extracted sources,
expanded all six generated function wrappers, and traced argument extraction,
return conversion, registration, and error handling.

- **RUSTSEC-2026-0176 (list/tuple iterator arithmetic):** Generated wrappers use
  FASTCALL and a separate borrowed tuple iterator for keyword parsing. Tuple
  arguments use fixed-index extraction; `read_ids` uses set iteration. DataFrame
  output conversion uses Rust column/vector iterators and constructs Python
  containers. The ordinary Python tuple iterator used in error formatting calls
  `next`, not the affected `nth` or `nth_back`. Our string `nth`/`skip` calls use
  Rust `str::Chars`, not PyO3 list/tuple iterators.
- **RUSTSEC-2026-0177 (missing `Sync` bound):** Function registration uses static
  `PyFunctionDef` values and `PyCFunction_NewEx`, not
  `PyCFunction::new_closure`. The reviewed Polars bridge and error conversions
  do not introduce that constructor. Internal Rust closures are not registered
  as Python-callable closures.

The inspected installed development extension also had no symbols for the
affected iterator methods or closure constructor/trampoline. This corroborated
the source review; symbol absence alone was not treated as proof.

#### Reassessment for nanalogue 0.2.1 (2026-10-23)

The [upgrade reassessment](https://ampcode.com/threads/T-01a168d6-aff0-75df-8eb8-e8c2db3785c5)
retained both exceptions for nanalogue 0.2.1 with the same Polars, pyo3-polars,
PyO3, ABI and extension-module versions/features listed above. The published
nanalogue 0.2.0 and 0.2.1 manifests have identical dependency requirements and
features. Its changed iterator paths operate on Rust slices, arrays and ranges,
not Python containers; it introduces no PyO3 calls. The binding's six public
signatures and their extraction, output conversion, registration and error paths
remain unchanged. Archive checksums and every extracted published regular file
were verified for nanalogue 0.2.1, PyO3, pyo3-macros-backend, pyo3-polars and
polars-error. The old nanalogue archive also matched its previous lock checksum.

On Rust 1.99.1, `RUSTC_BOOTSTRAP=pynanalogue cargo rustc --locked --features
pyo3/extension-module --lib -- -Zunpretty=expanded` confirmed all six wrappers
still use FASTCALL, with no tuple/dict wrapper or closure constructor/trampoline.
The rebuilt development extension installed for Python 3.11 matched the inspected
build's SHA256 and again lacked affected iterator/closure symbols; static
function registration and ordinary tuple `next` symbols were present.

The audit used cargo-audit 0.22.2 and cargo-deny 0.20.2 with RustSec database
revision `f7b2f6bc663452200bbe639ebb75e0c762b3422b` (2026-10-23). An isolated
lockfile scan with **no ignores** still reported only the two PyO3 advisories,
quick-xml's RUSTSEC-2026-0194/0195, and the bincode/custom_derive maintenance
warnings. `cargo tree --locked --all-features --target all -i` confirmed that
quick-xml and bincode remain inactive. **custom_derive is active** through
rust-htslib; its exception accepts an unmaintained dependency, not unreachable
code. With the existing exceptions, `cargo audit --json` reported no unignored
vulnerabilities or warnings, and
`cargo deny --locked check advisories bans sources` passed, also when explicitly
enabling `pyo3/extension-module`. No exception IDs were added or removed.
Upgrading nanalogue does not inherit its own published Cargo.lock or fix these
remaining advisories in this consumer's dependency graph.

This assessment covers PyO3 linked into this binding under the reviewed API and
features. It is not a formal whole-program proof or a security assessment of
Python Polars' separate native runtime or other extensions in the Python process.
ABI3 alone does not eliminate both vulnerabilities. Both advisories are fixed in
PyO3 >=0.29.0; the reviewed bridge and nanalogue Polars pin prevent a direct
compatible upgrade.

**Reassess these exceptions** when changing public argument/return types,
generated wrappers, callbacks, iterator adapters, dependency versions, or ABI/build
features. In particular, introducing Python list/tuple iterator skipping or
`PyCFunction::new_closure` invalidates the corresponding rationale and requires
review before retaining the exception. The tools match advisory IDs; they do not
automatically enforce this reviewed version/API scope. Remove the exceptions
when a patched PyO3 version is adopted.

## Disclosure Policy

Once a security issue is fixed:
- We will release a patch as soon as possible
- A security advisory will be published
- Credit will be given to the reporter (unless anonymity is requested)

## Contact

For security concerns, contact: mail AT unintegrable dot com

For general bugs and features, use the GitHub issue tracker.
