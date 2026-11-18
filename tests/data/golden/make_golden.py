"""Write the golden default outputs of ``scbamqc bam`` for the test fixtures.

The golden files must come from a version without the extra metrics, so that
``test_default_output_unchanged`` proves that the default output is the same.
They were written with the sources of commit ebfd952 (0.4.0.dev0):

    mkdir -p "$TMPDIR/ref" && git archive ebfd952 src | tar -x -C "$TMPDIR/ref"
    uv run --extra test --extra bam python tests/data/golden/make_golden.py \
        "$TMPDIR/ref/src" tests tests/data/golden
"""

import gzip
import os
import sys
import tempfile
from pathlib import Path

src, tests, out = (str(Path(p).resolve()) for p in sys.argv[1:4])
sys.path[:0] = [src, tests]

import scbamqc  # noqa: E402

assert scbamqc.__file__.startswith(src), scbamqc.__file__

from conftest import GOLDEN_RUNS, golden_fixtures, make_bam  # noqa: E402
from scbamqc.cli import app  # noqa: E402

tmp = Path(tempfile.mkdtemp())
os.chdir(tmp)
for name, (reads, contigs) in golden_fixtures().items():
    bam = make_bam(tmp / f"{name}.bam", reads, contigs=contigs)
    for run, args in GOLDEN_RUNS.items():
        cmd = ["bam", str(bam), "-j", "1", "--out-prefix", "o", "--no-progress"]
        try:
            app([*cmd, "--no-check-tags", *args])
        except SystemExit as e:
            assert e.code in (0, None)
        text = gzip.open(tmp / "o_per_cell.csv.gz", "rt").read()
        Path(out, f"{name}__{run}.csv").write_text(text)
