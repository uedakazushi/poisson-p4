# Actual N literal coneの独立レビュー保存

2026-11-04。独立算術と概念証明は
`AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63_REVIEW`の同時実行に保存した。
専用コードは`code/ambient_B18_N_same_R0_binding_review_20261104_v63.py`であり、
N source/R0結合とliteral coneを一度に監査する。このcompanionは同じ独立計算の
module-cone部分を公開する保存口で、重複計算やauthor再実行ではない。

八変数Koszulとv作用から36列・84列を独立に再導出し、全source relationと全鎖恒等式、
actual node代入後の全行列成分を確認した。正確性はR8上のKoszul解像度、
topologically nilpotentなTを使う単射v−Tとそのmonic coneから従う。
unit一個の消去後は最小で、Bettiは8,36,84,126,126,84,36,9,1。
実R0との結合、Cohen座標の許可、ν,a,zの入力は同時レビューで受理したものである。
D4 extension、T3、元B18原点の深さはこのレビューで解決していない。
