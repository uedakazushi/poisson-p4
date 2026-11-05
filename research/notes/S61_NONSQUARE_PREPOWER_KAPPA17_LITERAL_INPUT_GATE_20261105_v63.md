# κ17 の実入力切出しと全 phase guard 補足

2026-11-05。元 SOURCE_GATE の 75,390 owned は convolution/publication phase の上界であり、
元 checkpoint 全体を読む準備 phase の最大値ではない。この区別を旧 report の書換えなしで
本補足に固定する。旧 111,610 は matrix25x32 部分だけを数えた準備値だった。

新しい source-only 操作は、原三 JSON の係数 record を全て数え、20 atom の原 num/den を
一切変換せず別の小さな source 入力へ切り出すだけである。原三 JSON はそれぞれ
40、101604、1390 coefficient record を含む。切出し2858 record と metadata35,000を
同時に含めても140,892 slot以内。係数の和・積・多項式評価は行わない。

後続 worker は大きい checkpoint を再度読まず、この20 atom 入力、49個の既存 DAG 定義、
既存κ四多項式だけを読む。従って convolution phase の75,390上界を使える。
全原資料の準備 phase は本 source-only job に分離済みで、後続 worker の隠れた入力 cache
にはならない。source49の値とκ17は未計算である。

予測出力支持11728、整数 L1 高さ201 bit、denominator65 bitを保持する。
一時整数には2*(201+65)+16=548 bitを別に確保する。
新入力の整数・有理数・指数 record は元 JSON と文字通り一致させ、その source位置と
hashを保存する。scalar worker の許可は依然として別途必要である。
