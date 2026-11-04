# B7×3の35最大minorに対する費用preflight

入力は既にsource-boundかつgeneric rank3を証明済みの908項7×3だけ。行と列から
除く因子は既知a,d,F,B,C,Jと非零QQ定数だけに限り、全21成分で元値を再構成する。
さらにQQ column定数で整数係数に揃える。いずれもFULLU全域で可逆なので、
localized maxminor idealを変えない。

35個の三行組それぞれについて、六つのdeterminant termのmonomial支持集合を
Minkowski和で作り、そのunionを**支持上界**として保存する。係数の積和や
cancellationは一つも計算せず、新minorの零性・非零性を判定しない。次数と
整数係数bit長の保守的上界も入力だけから記録する。

このpacketはminor/GB jobではない。次の候補は、この小さくした整数入力の
固定35determinantsだけをFLINTで計算し、実支持数が十分小さい場合に別予算の
QQ saturation/Bezoutへ進む方法である。FULLU全域rank3やbadlocusの結論はまだ
出さず、fullH7がfullquartic239に必要であるとも仮定しない。
