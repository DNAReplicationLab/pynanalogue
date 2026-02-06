# Verifies that code examples in README.md produce their documented output.
# Uses HTML comment markers in README to identify testable code and expected output blocks.

import contextlib
import io
import re
from pathlib import Path

import pytest

README_PATH = Path(__file__).resolve().parent.parent / "README.md"

# Regex patterns for extracting marked blocks from README
# Matches: <!-- TEST CODE: START marker_name --> or <!-- TEST CODE: NOOUTPUT marker_name -->
_MARKER = r"([a-zA-Z0-9_]+)"
CODE_START_RE = re.compile(rf"<!--\s*TEST CODE:\s*(START|NOOUTPUT)\s+{_MARKER}\s*-->")
CODE_END_RE = re.compile(rf"<!--\s*TEST CODE:\s*END\s+{_MARKER}\s*-->")
OUTPUT_START_RE = re.compile(rf"<!--\s*TEST OUTPUT:\s*START\s+{_MARKER}\s*-->")
OUTPUT_END_RE = re.compile(rf"<!--\s*TEST OUTPUT:\s*END\s+{_MARKER}\s*-->")

# Matches markdown code fences like ```python, ```text, ```json, etc.
FENCE_RE = re.compile(r"^```\w*\s*$")


def _strip_fences(block: str) -> str:
    """Remove markdown code fences from a block, keeping only inner content."""
    lines = block.strip().splitlines()
    if lines and FENCE_RE.match(lines[0]):
        lines = lines[1:]
    if lines and lines[-1].strip() == "```":
        lines = lines[:-1]
    return "\n".join(lines)


def _normalize(text: str) -> str:
    """Strip trailing whitespace per line and trailing newlines."""
    lines = [line.rstrip() for line in text.splitlines()]
    return "\n".join(lines).strip()


def parse_readme_examples():
    """Parse README.md and extract all marked code/output pairs."""
    readme_text = README_PATH.read_text(encoding="utf-8")
    lines = readme_text.splitlines()

    code_blocks = {}  # marker -> (code_str, has_output)
    output_blocks = {}  # marker -> output_str

    i = 0
    while i < len(lines):
        line = lines[i]

        # Check for code block start
        code_match = CODE_START_RE.search(line)
        if code_match:
            block_type = code_match.group(1)  # START or NOOUTPUT
            marker = code_match.group(2)
            has_output = block_type == "START"
            i += 1
            block_lines = []
            while i < len(lines):
                end_match = CODE_END_RE.search(lines[i])
                if end_match and end_match.group(1) == marker:
                    break
                block_lines.append(lines[i])
                i += 1
            else:
                raise ValueError(
                    f"TEST CODE: START/NOOUTPUT '{marker}' has no matching END"
                )
            code_blocks[marker] = (_strip_fences("\n".join(block_lines)), has_output)
            i += 1
            continue

        # Check for output block start
        output_match = OUTPUT_START_RE.search(line)
        if output_match:
            marker = output_match.group(1)
            i += 1
            block_lines = []
            while i < len(lines):
                end_match = OUTPUT_END_RE.search(lines[i])
                if end_match and end_match.group(1) == marker:
                    break
                block_lines.append(lines[i])
                i += 1
            else:
                raise ValueError(f"TEST OUTPUT: START '{marker}' has no matching END")
            output_blocks[marker] = _strip_fences("\n".join(block_lines))
            i += 1
            continue

        i += 1

    # Validate: every code block with has_output=True must have an output block
    for marker, (_, has_output) in code_blocks.items():
        if has_output and marker not in output_blocks:
            raise ValueError(
                f"TEST CODE '{marker}' expects output but no "
                f"TEST OUTPUT block found for '{marker}'"
            )

    # Validate: every output block must have a code block
    for marker in output_blocks:
        if marker not in code_blocks:
            raise ValueError(f"TEST OUTPUT '{marker}' has no matching TEST CODE block")

    if not code_blocks:
        pytest.skip("No TEST CODE markers found in README.md", allow_module_level=True)

    # Build test cases: (marker, code, expected_output_or_None)
    cases = []
    for marker, (code, has_output) in code_blocks.items():
        expected = output_blocks.get(marker) if has_output else None
        cases.append((marker, code, expected))

    return cases


# Parse at module level so parametrize gets the IDs
_TEST_CASES = parse_readme_examples()


@pytest.mark.parametrize(
    "marker,code,expected_output",
    _TEST_CASES,
    ids=[c[0] for c in _TEST_CASES],
)
def test_readme_example(marker, code, expected_output, tmp_path, monkeypatch):
    """Execute a README code example and verify its output matches."""
    # Run NOOUTPUT examples (like simulate_mod_bam) in tmp_path
    # to avoid file pollution in the project directory
    if expected_output is None:
        monkeypatch.chdir(tmp_path)

    stdout_capture = io.StringIO()
    namespace = {}

    with contextlib.redirect_stdout(stdout_capture):
        exec(code, namespace)  # noqa: S102

    if expected_output is not None:
        actual = _normalize(stdout_capture.getvalue())
        expected = _normalize(expected_output)
        assert actual == expected, (
            f"README example '{marker}' output mismatch.\n"
            f"--- Expected ---\n{expected}\n"
            f"--- Actual ---\n{actual}"
        )
