#!/usr/bin/env python3
"""Run NFL_AB_design from a source checkout."""

from pathlib import Path
import sys


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PACKAGE_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from nfl_ab_design.workflow import main


if __name__ == "__main__":
    raise SystemExit(main())
