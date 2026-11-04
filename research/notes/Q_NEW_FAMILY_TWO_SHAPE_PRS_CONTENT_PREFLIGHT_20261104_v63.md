# 保存済みPRS prefixの既知contentだけを除くpreflight

2026-11-04。前の独立PRSはC1まで保存し、最終C0で2GB制限下の追加メモリ確保に失敗した。
この新しい一回60秒/2GB/32MiBジョブは**保存済みC2/C1**だけを読み、全係数に共通する
整数content、既知因子 `r,Z0,r−1,T=((r−1)R−2r)` だけを除く。
各完全除算と、最終 `C_i=F_i*primitive_C_i` の全係数再乗算を確認する。
未知の多項式gcd、PRSの再計算、最終norm、resultant、点、prime、source26は使わない。

削減後の次数・support・coefficient bit数と全配列を保存する。
このpreflightだけでは元の47,114係数raw resultant identityはまだ未検証である。

将来の別許可の最終段階は、primitive quadratic `a2*s²+a1*s+a0` と
primitive linear `b1*s+b0` に対する

```text
norm=a2*b0²−a1*b0*b1+a0*b1²,
C0=(F1²/F2)*norm/a2².
```

である。これはここでは評価しない。必要なメモリの見積もりは保存したsupportから行い、
全同一性や追加計算の成功を先取りしない。
