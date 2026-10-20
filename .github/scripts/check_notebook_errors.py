# ABOUTME: Fail CI when an executed notebook contains an error, including tracebacks
# ABOUTME: rendered as display/stream output (rich handlers, `!cmd` CLIs) not just `error`

"""Usage: python check_notebook_errors.py <executed.ipynb>

Exits 1 if any cell output is an ``error`` or contains a Python traceback,
printing the offending cells. Tracebacks are searched for in stream text and in
``text/plain`` / ``text/html`` display data after stripping ANSI codes and HTML
tags, since rich prints ``Traceback`` and ``(most recent call last)`` in
separately styled spans.
"""

import json
import re
import sys

ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
TAGS = re.compile(r"<[^>]+>")
TRACEBACK = re.compile(r"Traceback\s*\(most recent call last\)")


def _texts(output):
    if output.get("output_type") == "stream":
        yield output.get("text", "")
    for key in ("text/plain", "text/html"):
        value = output.get("data", {}).get(key)
        if value is not None:
            yield value


def _as_str(value):
    return "".join(value) if isinstance(value, list) else value


def main(path):
    with open(path, encoding="utf-8") as f:
        nb = json.load(f)

    bad = []
    for number, cell in enumerate(nb.get("cells", [])):
        for output in cell.get("outputs", []):
            if output.get("output_type") == "error":
                bad.append((number, f"{output.get('ename')}: {output.get('evalue')}"))
                break
            text = " ".join(
                TAGS.sub("", ANSI.sub("", _as_str(t))) for t in _texts(output)
            )
            if TRACEBACK.search(text):
                bad.append((number, "traceback in cell output"))
                break

    for number, reason in bad:
        source = _as_str(nb["cells"][number].get("source", "")).strip().splitlines()
        first = source[0] if source else ""
        print(f"{path}: cell {number}: {reason}\n    {first}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
