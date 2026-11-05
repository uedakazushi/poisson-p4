# Exact43 mixed-source binary compaction

The current43 residual is193 by40 with400727 logical terms. The source plan resolves every active JSON row or binary snapshot cell by explicit encoding, offset, length and SHA256. Each source value is independently encoded/decoded; no full coefficient matrix is resident. The output snapshot uses the accepted uncompressed exact signed-integer/exponent schema and retains all original row/column identifiers including zero rows.

Fresh actual-read inputs plus every new plan/snapshot/index/code/IO/equality record must fit32MiB RAW, owned150000 slots,60s/2GB. All earlier text/binary generations and stop records remain stored. No coefficient arithmetic, pivot, nonmonomial inverse, GCD, numeric point, prime, rank or exhaustion is part of this representation stage. The gate is359fe767a2712ce986754e29366d56d4ff14c3b2e53b5947356fe5093a710556; compaction must verify every source/output value and row roundtrip before atomic publication.
