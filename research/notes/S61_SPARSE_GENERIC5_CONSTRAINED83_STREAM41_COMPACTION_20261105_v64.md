# Exact current41 compaction

Under source gate7e16b0e86e5ae321aed2ec7fdd2c00f3f05c1348c97bacf174f66d39468c665f this copies the same195 by42 current41 matrix,301097 terms and4853 cells. Every literal payload is source-hash checked, copied byte for byte, read back twice and compared; all row identifiers including zeros persist. The fresh actual-read ledger and all new exact snapshot/index/receipt/code data fit RAW32MiB, owned150000,60s and2GB. Every older generation remains stored. No polynomial arithmetic, pivot or new inverse is performed.
