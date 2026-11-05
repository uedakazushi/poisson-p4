# 列17の pair 補正第四成分に対する source-only gate

2026-11-05。corrected 四＋十分解の最初の実 scalar 候補を一つだけ固定する。
A4逆の最後の行は κ/T なので、列17の四 pair 成分 v17 に対し

    x4 = (κ0 v17,0 + κ1 v17,1 + κ2 v17,2 + κ3 v17,3) / T.

κ の四多項式は既存の独立 QQ 確認済み root179/198/217/237 をそのまま用いる。
符号をもう一度掛けない。v17 は current15 [2,4,5,8] の actual column17。
この列の normalization は1で、Q2/L source DAG の49 node/20 atomを文字通り使う。
全 denominator は QQ 定数と a の単項式で、B の払いは不要である。
κ17=κ dot v17 の分子だけを予定し、他列、他の三つの adj4 行、最終 r17 を計算しない。

source support は既存 atom の box/全次数/QQ L1 上界を和・積で伝播する。
κ の既存値から同じ上界を作り、四積を一つずつ生成して accumulator に加える。
49 source 値は last-use で解放し、四つの v17 を保持する。
全 resident input、source cache、積、更新前後 accumulator、publication copy、
metadata と整数 temporary を含む上界を report に保存する。
これは SOURCE gate のみで、実 polynomial worker は親の別の許可を要する。

κ17 が非零でも、全1×8や cubic57 を証明しない。零なら列17について x4=0 とだけ
言え、raw pair 比例や他の補正の零化は従わない。四＋十の残りの補正を全て保持する。
