# B18 の五有理消去・十二 flag 商式の独立監査

PROPER_RATIONAL_ELIMINATION の凍結 packet を、作者の builder を使わずに別実装の
sparse QQ 多項式辞書と分母の交差積だけで確認する。symbolic cancel、GB、新 source
jet、数値点、次段の消去は実行しない。

元の14変数17式に対し、順序 F2(v0), F0(b), F4(v1), F3(e), F1(c) の五恒等式を
全係数で検査する。各先頭係数は 1/3,2,2,-(1-d²),2。D=1-d² の原点値は1なので、
selected5 の ideal は明示された五有理根の graph ideal と局所環・完成局所環で等しい。
単に根を代入した零性だけから ideal の等号を推測しない。

残り十二の元 Fi について、五有理根を独立代入し、全 rational numerator と constant
denominator を確認する。保存された整数 numerator が literal に content*d*Hi と等しい
ことを全係数で検証する。従って、先の独立受理済み source-gauge 模型を通じ、全 E4 は

    QQ[[a,d,v2,...,v8]]/(d*H0,...,d*H11)

である。原点12×9 Jacobian と指定5×5 minor 19660800も商多項式から別計算する。

別途 divided ideal (H0,...,H11) を考えた core と、全 E4 の ideal (dHi) を混同しない。
Hi の selected5 を全 E4 で零にしてはいけない。actual D0、CM、depth、socle、次段
有理消去は本監査の結論ではない。結果・時間・入力 hash は report に保存する。
