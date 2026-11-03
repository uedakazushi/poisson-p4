# 9×43 criterion の成分数表記の訂正

2026-11-03。凍結済み `S61_NONSQUARE_CUBIC_NINE_ROW_CRITERION_20261103_v63` の note と
引継ぎ本文で「全5190成分」と記した箇所は、**57×91=5187個の行列成分**が正しい。
matrix file の header は3整数だが、それらは行列成分ではない。
単なる積の表記誤りであり、native/Python の実際の loop は最初から57行×91列で、
`len(lines)==1+57*91` を検査していた。

行列 shape57×91、係数項数164333、heavy/pure分割、9×43の387成分、49角 bordered minor、
全 source hash、点別 QQ 検査と数学的結論は変わらない。旧 packet の bytes は変更しない。
