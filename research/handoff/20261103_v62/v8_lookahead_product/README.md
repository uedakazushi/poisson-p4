# Frozen V8 lookahead checkpoint

The exact scope and mathematical argument are in
`research/notes/S61_BINARY_NONSQUARE_LOOKAHEAD_20261103_v62.md`.
This packet does not prove all-point product/cubic rank, quartic65, or ambient exhaustion.

Extract the archive at repository root. Reuse the published Cloud activation:

```
source /workspace/cloud-setup/activate.sh
w=work/nonsquare-v62/product-field-v7
c++ -O2 -std=c++17 $w/low-image-lift.cpp -o /tmp/s61-v8-low-image $(pkg-config --cflags --libs flint)
/tmp/s61-v8-low-image $w/row1-direct-module-source 3 /tmp/s61-v8-k10.kernel $w/row1-direct-D3.bounds 10
cmp /tmp/s61-v8-k10.kernel $w/low-image-D3-K10.kernel
```

The saved `witness-image-stats.cpp` replays the complete image, verifies a known-unit target,
adds/replaces its column, and performs exactly one native unit pivot. Its CLI is
`input output units metadata01 witness targetrow`. The source input and six units are
pinned public PRODUCT_ALT archive members; hashes appear in the manifest. For k10 use
`low-image-D3-K10.witness`, targetrow1, metadata1. Expect `EXACT_WITNESS_APPEND 46`,
13×46 and167410 terms. For k12 use its separate witness; expect a true replacement,
13×45 and186299 terms. Do not merge these two presentations.

Large generated matrices are omitted with SHA256 and exact replay instructions. The main
full-cubic27×61 matrix remains in immutable V6. This packet's thin product case is separate.
The independent BOUNDED_UNIT_LIFT review predates V8; it is not a review of the lookahead.
No solver owned by this agent was running at the Cloud outage or restoration audit.
