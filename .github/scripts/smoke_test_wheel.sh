# Install one existing release wheel in isolation and run its functional tests.
# No Rust toolchain, local package build, or publication is involved.
set -euo pipefail

wheel_dir="$1"
source_dir="$2"
shopt -s nullglob
wheels=("$wheel_dir"/*.whl)
if [ "${#wheels[@]}" -ne 1 ]; then
  echo "Expected exactly one wheel in $wheel_dir, found ${#wheels[@]}" >&2
  exit 1
fi

scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT
"${PYTHON:-python}" -m venv "$scratch/venv"
python="$scratch/venv/bin/python"
export PYTHONDONTWRITEBYTECODE=1
# Refuse source builds, including for runtime dependencies.
"$python" -m pip install --only-binary=:all: "${wheels[0]}" pytest
"$python" -m pip check
"$python" - "${wheels[0]}" "$source_dir/pyproject.toml" <<'PY'
import hashlib
import importlib.metadata
import platform
import sys
import tomllib
from pathlib import Path

import pynanalogue

wheel = Path(sys.argv[1])
expected = tomllib.loads(Path(sys.argv[2]).read_text())["project"]["version"]
actual = importlib.metadata.version("pynanalogue")
assert actual == expected, (actual, expected)
print(f"Wheel: {wheel.name}\nSHA256: {hashlib.sha256(wheel.read_bytes()).hexdigest()}")
print(f"Python: {platform.python_version()} / {platform.machine()} / {platform.libc_ver()}")
print(f"Installed pynanalogue {actual}: {pynanalogue.__file__}")
PY

cd "$source_dir"
"$python" -m pytest tests/ --ignore=tests/test_benchmark.py -v \
  -p no:cacheprovider --basetemp="$scratch/pytest"
