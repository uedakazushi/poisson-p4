# 重み6閉切断の二軸idealの独立監査

2026年11月3日。親37模型の全有限DAGを別のFraction辞書実装で代入し、
生150係数と全350残余の全係数を再構築した。保存した座標変換を行い、
350×6のforwardおよび6×350のreverse恒等式をQQ[a,W1,W2,eta]で全て掛け戻した。
reverseはa倍の六生成元を与えるため、必要な逆元はaだけである。
Smith routineの対角出力に全ideal等号を依存させていない。

独立な小さいM2計算で、三つのprimary idealの交叉、nilradicalの平方零性、全次元1も照合した。
主分解がminimal三成分だけを持つのでtotal閉切断にはembedded primeがなく、
一次元の全total schemeはCMである。nilradicalはa軸上の自由3加群で、被約支持は
二つのspecial値で交わる曲線。各crossingの完成環は

    QQ[[u,W1,W2,eta]]/(W1²,W1W2,W2²,W1 eta,W2 eta,u eta²)

であり、depth1・multiplicity5を小M2で確認した。
一方、そのspecial fiberは(W1,W2)*(W1,W2,eta)をidealとし、depth0。
totalにembedded primeがないことをspecial fiberへ誤って拡張しない。

既存T40模型との負軸全150係数一致は、入力不変のfresh producerで再現した。
正軸の全成分帰属、全親R37の新しい非平坦性、全37変数の枝の尽くしは承認していない。
保存先は `code/q_weight6_two_axes_review_20261103_v63.py` と
`reports/Q_GLOBAL_WEIGHT6_TWO_AXES_20261103_v63_REVIEW/`。
