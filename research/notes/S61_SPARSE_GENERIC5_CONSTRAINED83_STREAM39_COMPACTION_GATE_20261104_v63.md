# Source-only exact-current compaction gate

The current 39-pair generation is 197 by 44 with 234898 logical polynomial terms. This gate resolves each row reference into its actual literal files and verifies the active immutable row receipts. It indexes one cell at a time, pins all source payloads by SHA256, and enumerates every file that the fresh compaction job will actually read, including current metadata, index and receipts.

The existing 23923327-byte retained-history accounting remains intact. Those retired generations stay stored unchanged. The parent's new fresh-job budget counts the explicit actual-read input set plus the new complete literal snapshot, new index and all new code, IO records, metadata and receipts. This distinguishes inactive retained history from actual current input; it performs no deletion and makes no claim that old bytes disappeared.

Compaction, if this gate passes, copies every source polynomial payload byte for byte into one complete sparse matrix file, independently reads it back, and verifies source/output byte equality and hashes for every cell and row. All zero rows and original identifiers remain. A new index references exact output offsets. Its owned-memory proof includes current/new metadata, a bounded cell buffer and cell-parse scratch. Logical full-matrix terms remain external. No polynomial arithmetic, normalization, row operation, GCD, point, prime or new inverse is required or authorized by compaction.

The gate uses the unchanged limits of 60 seconds, 2 GB, 150000 owned scalar slots and 32 MiB RAW for the fresh actual-read input plus all generated artifacts. It does not replay any old pivot or assert another unit pivot, a completed left inverse, or exhaustion.
