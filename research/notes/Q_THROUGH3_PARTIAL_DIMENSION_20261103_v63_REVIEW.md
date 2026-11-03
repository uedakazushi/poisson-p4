# 小さい13式証明書と全R35の上界22：独立監査

2026年11月3日。全74式の新たなGBは実行していない。
先に独立監査済みのQQ57二次式・17真三次式から、32003剰余体の入力を全係数で再構成し、
保存 input.m2 と文字通り一致することを確認した。全分母はpで単元である。
全report/source/evidenceと解凍後hashも確認した。

保存された小さい13列の実 cofactor 恒等式 A*C=B を新しいM2プロセスで検査し、
**実際の**13 leading monomials が保存一覧と一致することも検査した。
全329列の追加再生は今回は行っていない。元の別プロセス22秒再生の凍結証拠と実装を監査した。

13先頭単項式のsupportは非空で二つずつ素である。
従って product criterion から選択した13多項式はそれ自体でGBをなし、そのidealは高さ13、次元22。
その全てがspecial J74に属するという実 polynomial identity を確認しているため、special J74 の次元は22以下。
また全329保存 leading supports を高さ13のcoordinate primeが覆うことも確認した。
これと13本の完全交差を合わせ、その部分leading idealの次元は正確に22となるが、
全 initial ideal がこの部分idealと等しいとは言っていない。

標数零への移送も受理する。実際のQQ homogeneous74式をV=Z_(32003)上で定義し、
Proj(V[u,z]/J74)のpropernessによるfiber dimension上半連続性を使う。
special affine cone≤22、したがってspecial Proj≤21からgeneric Proj≤21、generic cone≤22が従う。
空Projの場合もcone dimension0で同じ上界となる。
全初期化と剰余体base changeの交換、全J74のflatnessは仮定していない。

J74_QQは真のm-adic initial idealに含まれるという先のQQ証明から、
全R35の associated graded と完成環の次元の等しさにより、複素数上の全完成局所モデルもdim≤22。
これは単なる有限体サンプルのdimをgenericdimへ移した議論ではない。
目標12・正確なdim22・全枝分解は未証明であり、degree4の3GiB制限停止もそのまま保持する。

再現コードは code/q_through3_dimension_review_20261103_v63.py、
小さいM2 traceと結果は reports/Q_THROUGH3_PARTIAL_DIMENSION_20261103_v63_REVIEW/ に保存する。
