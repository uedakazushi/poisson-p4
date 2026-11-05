# Exact47 active-fragment mixed binary compaction

The current47 matrix is unchanged. The old8,028,307-byte input snapshot is still stored, but its265 cells still referenced by current47 are independently byte/value copied into a52,877-byte fragment, with complete original-index coverage and hash map87f60dce... . Current new binary row/index pairs remain unmodified. This source plan covers every4094 current cell and189 rows using the actual fragment plus all new rows; the old superset file is not read by this fresh compaction.

Every old byte stays stored; nothing is deleted or hidden. Each actual-read file, map, new snapshot/index/code/IO/receipt is freshly RAW charged. Hard32MiB/150000-owned incl overlay metadata/60s/2GB unchanged. No pivot, coefficient arithmetic, GCD, point, prime, nonmonomial inverse, rank or exhaustion occurs here.
