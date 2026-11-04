# A₂−kA₁²の最高Fn極：定数kによる相殺はない

2026-11-04。二shape写像のA₁=N₁/(50Fₙ)、A₂=N₂/(400Fₙ²)について、**定数k∈QQでA₂−kA₁²の分子をFₙで割り切れるようにはできない**ことを、小さい完全係数検査で確認した。cubic/quinticの二式はこの案では小さくならなかった。

実際にP=N₂/400、Q=N₁²/2500をQQ[r,s,R]の多項式として扱い、単一の既知因子Fₙで割った余りをそれぞれ求めた。各除算は商×Fₙ+余りが元に戻ることを確認した。二つの余りの係数vectorは比例せず、二つの係数行における非零2×2 minorをcertificateへ保存した。従ってP−kQがFₙの倍数となるQQ定数kは存在しない。

これは単一多項式による既知因子の完全除算・係数比較であり、generic multivariate gcd、GB、resultant、source26、新しい点を使わない。kをQQ(R)、QQ(r,R)の有理関数に拡大した案や、別のtarget多項式の組合せまで排除するものではない。自動的に大きなansatzへ進まなかった。

producerは`code/q_new_family_two_shape_pole_reduction_20261104_v63.py`、保存は
`reports/Q_NEW_FAMILY_TWO_SHAPE_POLE_REDUCTION_20261104_v63/`。
規格化済みの全分子、全余り、候補定数の不一致、二係数minor、入力hashを保持する。計算約0.087秒。前の二shape式のbidegree(9,3)、(14,5)、generic有理有限性、未知target体index・方程式の状態は変わらない。
