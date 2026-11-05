# Q1 の J 付値0：独立 ROOT d-Horner 検証

2026-11-05。独立報告
`S61_NONSQUARE_PREPOWER_Q1_J_ROOT_DHORNER_REVIEW_20261105_v63`
の SHA は `d54591971ecf7db1f6670f40d1783ac99f0b591e05b2376268ae3a986f21cfc4`。

ROOTは元の4811整数をd次で読み直し、通常の下降d-Horner法でk0の整数層を再構成した。
producer imports、binomial係数、root-power cacheを使わず、全606係数、numeric D5¹⁸と
symbolic a¹⁸が一次結果と一致した。19213 products、343 bit、28046 owned slots。
これにより普通の \(\mathbf Q[a,e,d]\) で \(v_J(Q_1)=0\) を受理する。

一次報告・準備資料は変更しない。J商、他因子やGCD、二tail idealの互いに素、
十分開の補集合消滅、全FULLU57/239は主張しない。
