# ABOUTME: Importing jaff must not replace the interpreter's exception hook; a global
# ABOUTME: rich traceback handler hid notebook errors from CI (display_data, not error)

import subprocess
import sys


def test_import_keeps_sys_excepthook():
    """``import jaff`` leaves ``sys.excepthook`` untouched.

    Runs in a fresh interpreter so no earlier import in the test session can
    mask the side effect.
    """
    code = (
        "import sys\n"
        "before = sys.excepthook\n"
        "import jaff\n"
        "assert sys.excepthook is before, sys.excepthook\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
