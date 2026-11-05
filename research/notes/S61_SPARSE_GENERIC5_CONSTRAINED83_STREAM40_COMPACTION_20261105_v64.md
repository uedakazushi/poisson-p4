# Exact current40 compaction

This copies the accepted current40 residual,196 by43 and292652 logical terms, into one literal text snapshot and plain compact-JSON index. Gateccca04eef129b9a99a6ef14088f8d6016a770fdcd52d545d08660b4f06617d20 proves the fresh actual-read input plus every new snapshot/index/receipt/code artifact fits RAW32MiB. Current and output metadata plus one cell remain below150000 owned scalar slots; time60s/memory2GB.

Each of4949 source payloads is hash-checked, copied literally and independently read back. Every196 row identifier including zero rows is retained, with per-row source/output equality receipts. No polynomial arithmetic, pivot or inverse is performed. All previous files remain stored unchanged; the fresh ledger explicitly lists actual-read inputs. Smaller JSON formatting is uncompressed exact metadata, with no coefficient loss or cap increase.
