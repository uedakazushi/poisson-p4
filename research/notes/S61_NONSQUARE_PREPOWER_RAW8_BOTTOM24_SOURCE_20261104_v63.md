# 零raw8の同じ24 source係数だけの抽出

固定raw8 workerは全15/90寄与まで完成し、選択minorは零となった。node auditでは
六つのheavy2 minor中五つが非零、選択4×6pure部分の全90二次minorが零である。
この六列以外のpure rankや全pair-kernelについては結論しない。

独立な零証明を小さくするため、同じ凍結native-inputから24 source entryだけの190node/
46atom closureを抽出する。1898 workerもdeterminantも再実行しない。全190値を保持する
過大見積81760、scratch5*640とreserve15000を加えて99960項以下。32MiB/2GB/60秒を保つ。
原子は既存inputの逐語recordを使い、新しいsource lookupやcoefficient queryをしない。

保存する4×6はcurrent15 rows `[2,4,5,8]`、pure columns元32
`[22,11,13,25,15,19]`。この四行のa正規化冪は全て0である。
独立peerはこの24polynomialと元source bindingを照合し、全90個の2×2積差を確認できる。
全17pure columnsへrank一を外挿しない。特に全8Q9tailを捨ててよいことにはならない。
