"""Build a platform wheel bundling the msentity-similarity C library."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
LIB_DIR = HERE / "msentity_similarity" / "lib"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", required=True, type=Path)
    parser.add_argument("--platform-tag", required=True)
    parser.add_argument("--out", type=Path, default=HERE / "dist")
    args = parser.parse_args()

    # Read through a versioned .so symlink, but retain the public filename
    # (libmsentity_similarity_c.so) that the runtime loader searches for.
    native_name = args.native.name
    native = args.native.resolve()
    if not native.is_file():
        parser.error(f"native library not found: {native}")
    args.out.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(LIB_DIR, ignore_errors=True)
    LIB_DIR.mkdir(parents=True)
    shutil.copy2(native, LIB_DIR / native_name)
    try:
        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run(
                [sys.executable, "-m", "pip", "wheel", str(HERE), "--no-deps",
                 "--no-build-isolation", "-w", tmp], check=True)
            (built,) = Path(tmp).glob("msentity_similarity-*.whl")
            subprocess.run(
                [sys.executable, "-m", "wheel", "tags", "--remove",
                 "--python-tag", "py3", "--abi-tag", "none",
                 "--platform-tag", args.platform_tag, str(built)], check=True)
            for wheel in Path(tmp).glob("msentity_similarity-*.whl"):
                shutil.move(str(wheel), args.out / wheel.name)
                print(args.out / wheel.name)
    finally:
        shutil.rmtree(LIB_DIR, ignore_errors=True)
        shutil.rmtree(HERE / "build", ignore_errors=True)
        for info in HERE.glob("*.egg-info"):
            shutil.rmtree(info, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
