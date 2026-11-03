# Frozen V9 short unit pivots

See `research/notes/S61_BINARY_NONSQUARE_SHORT_UNIT_PIVOTS_20261103_v62.md` and STATUS.json.
Both cases are invertible column replacements, followed by one known-unit pivot.
There is no new all-point rank, quartic65 or ambient conclusion.

Extract the archive at repository root. Reuse `source /workspace/cloud-setup/activate.sh`.
The generic search source `primitive-unit-span-optimized.cpp` has CLI:

    input known-units metadata01 target-row output-witness max-unit-factor-count

Compile with `c++ -O2 -std=c++17` and `pkg-config --cflags --libs flint`. The final main
witness can be rediscovered from the included one-row61-column main-target-row.txt,
using P-units, metadata0, targetrow0 and factor count8. It exactly matches the saved
61-component witness. The product case can use its full14×46 input, metadata1,
targetrow10 and factor count8 (the first version already found it at count3).

The included witness-image-stats.cpp has CLI:

    input output known-units metadata01 witness target-row

It prints the complete image statistics, verifies the target unit, makes the column
replacement and performs exactly one native pivot. Product expects replacement11,
13×45/44143 terms. Main expects replacement5,26×60/1290701 terms; reconstruct P-clean27
from immutable V6 first and verify its SHA from this manifest. In the main serialized
three-variable format the third variable named d represents P=ad.

Large generated outputs are omitted with SHA256. The two original target rows and
full small product source are included so the unit identities can be independently
checked without rebuilding large matrices. No solver owned by this agent was left
running after these replays.
