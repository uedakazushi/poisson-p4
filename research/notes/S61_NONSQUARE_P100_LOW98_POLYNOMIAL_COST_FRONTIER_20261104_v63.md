# P100 low98 polynomial matrix の cost frontier

同じ98 sourceを使う実 global minorを固定できたが、その coefficient/rank 演算を実行する前の resource gate が残る。保存済み PURE70 forecast は70行を同時保持すると454974 RF terms、最も大きい一行の streaming bound は56407。七つの98対84 multiplication operatorを全て実体化するだけでも forecast134032 termsとなる。この三つの数字は既存 inputから再読して照合した。

150000 live termsの上限内で一行ずつ coefficientを生成することはできる。しかし polynomial rank basis、determinant、global gcd の stateを同じ上限へ収める証明はまだない。従って resident98行列や polynomial Gaussian は開始しない。この記録は全ての低memory algorithmが不可能という結論ではない。

一つの有効な将来の positive criterion は、実 low-source determinant の numerator Dを known-unit denominator付きで計算し、Dのmonic化が32003-integralなunit leading coefficientを持つと確かめ、その有限 quotientの全剰余体で実sourceのrankを閉じる方法である。H15で成功した有限CRT/Nakayamaと同じ論理だが、新Dの値・次数・factor degreeは未計算である。H15の七つの成功をそのままDの全零点へ外挿しない。

新しい normal form、coefficient、minor、rank、点、prime、GBは実行していない。全P100239は未解決。
