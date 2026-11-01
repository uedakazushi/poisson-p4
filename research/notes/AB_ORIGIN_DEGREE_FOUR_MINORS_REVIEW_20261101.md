# AB 四次：全249評価minorの別実装検算

2026-11-01。対象は `AB_ORIGIN_DEGREE_FOUR_20261101` の下界38570である。

主producerが保存した45座標多項式、726個の十変数評価点、249個のweight群、
各群で選んだ四次積から、**全ての正方評価行列を別に再構成**した。
Pythonの `Fraction` で有理係数を還元し、NumPyの整数演算で座標と積を評価する。
各座標の全項が指定weightを持つこと、各四次積が群のweightを持つこと、
同じ群の選択積が重複しないことも検査する。

正方行列には合計15,561,408成分がある。
各乗算の直後に32003で還元するため、NumPy int64の乗算にはオーバーフローがない。
別native verifierはFLINTを使用せず、行交換とFermat逆元による通常のGaussian消去で
行列式を計算する。32bit積の上限32002²は2³²未満である。
全249行列式が主producerの非零値と一致し、独立列数の総和は38570となった。
構成した全行列のbinary streamのSHA-256もreportに保存した。

fresh実行は1.401秒でPASSした。

```sh
.venv/bin/python code/ab_origin_degree_four_minors_review_20261101.py --output /tmp/ab-fourth-minors-review.json
```

これは保存された座標多項式からの独立評価・全minor検算である。
incidenceのQへの写像自体を第二の幾何学的構成で導いたものではない。
cohomologyによる上界38570、P-cyclic評価、全次数正規性の証明とも区別する。
異なるweightの単項式支持が互いに素であることと、全係数分母が32003で可逆であることにより、
この証明書はQQ上の座標四次の独立性を与える。
