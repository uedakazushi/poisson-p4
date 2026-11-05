# Exact active subset of retained binary source

The current47 manifest uses newly computed binary row/index pairs plus only a subset of its former input snapshot. This representation stage selects exactly every old-source cell index still referenced by that manifest, copies its original binary payload byte for byte with independent decode/re-encode and readback checks, and publishes a complete replacement map. Current new-row files are neither touched nor omitted.

All old snapshot bytes and every earlier generation remain stored unchanged. A later whole-current compaction may read the smaller exact fragment instead of the old superset; its ledger must include the fragment, replacement map, current new rows and all metadata. This is not deletion, coefficient loss, precision change, gzip or a new matrix operation. Each fresh job is separately preflighted under32MiB RAW,150000 owned incl metadata,60s/2GB. No pivot, GCD, point, prime, nonmonomial inverse, rank or exhaustion is involved.
