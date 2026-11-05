# graph5 V1 の値実行前 SOURCE-HOLD 記録

2026-11-05。独立precheck report
`reports/AMBIENT_B18_C3_GRAPH5_CACHED_SOURCE_PRECHECK_20261105_v63_REVIEW/report.json`
SHA256 `a56c4a4d0d848a132ebb9357b4d75a83f4850620002a3580533a5a7614db0816`
は `all_passed=false`、値実行の許可条件は未充足である。

凍結された元graph5 SOURCE/VALUE V1は、op-cの分母r[3]を読まずr[2]だけを使っていた。
45定数中32個が非整数で、例としてnode20084の1/4が1、node1675の−4/5が−4になる。
compiled DAGそのものは元sourceを保っている。
また五列のbA2には実際にmixed H項があり、「全M2がzzだけ」という文は誤りだった。
eB射影の数学的な修復にはr/ℓ倍RHによるclearingと、RH2_eB-plane=0の使用が必要である。

全32定数と五行の実mixed項を別証拠に保存した。rootや三次係数のreadbackは一切していない。
元V1は書き換えず、別V2に修正を残した。後続のFORCED_ZERO純粋証明により、
五係数の消滅は数値実行なしに証明された。
