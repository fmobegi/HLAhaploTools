"""Set the package version in pyproject.toml (used by the TestPyPI workflow for its .devN builds).

python scripts/set_version.py 0.1.0.dev3
"""

import re
import sys
from pathlib import Path

version = sys.argv[1]
if not re.fullmatch(r"\d+(\.\d+)*((a|b|rc)\d+)?(\.post\d+)?(\.dev\d+)?", version):
    sys.exit(f"error: {version!r} is not a PEP 440 version")
path = Path(__file__).resolve().parent.parent / "pyproject.toml"
text = path.read_text(encoding="utf-8")
new, n = re.subn(r'(?m)^version = "[^"]*"$', f'version = "{version}"', text, count=1)
if n != 1:
    sys.exit("error: no version line in pyproject.toml")
path.write_text(new, encoding="utf-8", newline="\n")
print(f"version {version}")
