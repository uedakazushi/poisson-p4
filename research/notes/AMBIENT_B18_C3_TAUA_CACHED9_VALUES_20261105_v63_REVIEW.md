# 九 tauA cached cubic 射影の独立 VALUE 検算

2026-11-05。SOURCE peer `b7cd1790…` 通過後の一回の固定 cache job を検算した。
PRIMARY report は `AMBIENT_B18_C3_TAUA_CACHED9_VALUES_20261105_v63/report.json`
SHA256 `0bc5fcb9d6e4343a4436192a298d79cd2843e82281564c7da3637ab0acfb65a3`。
独立 REVIEW report は
`reports/AMBIENT_B18_C3_TAUA_CACHED9_VALUES_20261105_v63_REVIEW/report.json`
SHA256 `fb749f09987ee2c9bc62c3a40f8346dc6b285fb135d6f19d70df852a22ea5d47`。

sourceの12×9 kernelがliteralなfree-column embeddingであることを全108項で再確認した。
独立 parser で元C2のdelta係数index29・30を読み、順序を逆にして

\[
-\frac8{25}\delta_{r^3}+\delta_{r^2\ell}
\]

を計算した。全九値は正確に0であり、PRIMARYの保存した18個のterm productも全て一致した。
選択外の三pivot行の追加pairingや、source/gauge/物理jetの再計算は行っていない。
Fractionの積・和の前にbit上界を確認し、処理は約0.0018秒で完了した。

これは九個のeA-r²ℓ射影に限るnegative結果である。
他の三次方向・quartic補正は残り、C3上界19、長さ18、torsion指数の改善を導いていない。
authorコード・SOURCE・VALUEのfrozen bytesは保持した。
