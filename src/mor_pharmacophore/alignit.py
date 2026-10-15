"""Step 3: Align-it scoring (external binary, GPL-3.0, https://github.com/OliverBScott/align-it)."""
import shutil
import subprocess


def alignit_path():
    return shutil.which("align-it")


def alignit_help():
    exe = alignit_path()
    if exe is None:
        raise FileNotFoundError("align-it is not on PATH")
    return subprocess.run([exe, "--help"], capture_output=True, text=True)
