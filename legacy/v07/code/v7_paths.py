"""Portable paths; expected certificates remain unchanged on reproduction."""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT/'data'
DATA=Path(os.environ.get('P4_V7_OUTPUT_DIR',str(ROOT/'data'/'reproduced'))).resolve()
DATA.mkdir(parents=True,exist_ok=True)
