# Exact current40 compaction source gate

This source-only gate resolves the accepted exact40 manifest (196 by43, 292652 logical terms) to its actual active literal files, current index and row receipts. All previous generations remain stored unchanged. The fresh-job ledger counts every file actually read and all prospective snapshot/index/receipt/code output, under60 seconds,2GB,150000 owned scalar slots and RAW32MiB.

A current-cell index is written as plain compact JSON, without compression. At most one source polynomial is parsed at a time. The forecast covers an exact full text snapshot, exact metadata index and a1MiB allowance for all additional source-gate/worker/IO/equality records. No polynomial update, pivot, GCD, point, prime or old-step replay occurs. A failed forecast does not mean all representation methods or units are unsafe.
