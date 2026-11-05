# 正規化済みC四乗商のF-Hasse：未実行worker

2026-11-05。ROOT integer longdivision peerの出力をそのまま再利用する。
ここでQ_C=T1/(a³C⁴)であり、原S61のquadric q1ではない。
独立検証済み4811整数係数／197-bit Dのbodyをbyte同一でコピーし、headerだけ変える。
scalar正規化、商、F power、Hasse層は新たに計算していない。compileとhelpのみである。
fresh ROOT/source-code監査と一回GOを実行前に要する。

F=500a+90e²+3e−24、a0=(24−3e−90e²)/500。
`R=24−3e−90e²`のl1は117<500、入力integer l1は2³⁰¹未満。
a次数16なのでcommon D*500¹⁶の整数Hasse係数をk=0..16でstreamし、最初の非零で止める。
最初の非零層だけ逆順に引いて零へ戻し、再現して保存する。
全integer partialはinput_l1*2¹⁶*500¹⁶未満、463 bits以内で、hard518は上げない。

入力4811、root cache289、e≤65,d≤18のaccumulator1254を持つ。
coefとpacked key両方およびmetadata/scratch16000によりwhole28708 cells。
input vectorは4811を先にreserveし、source L1もcommit前transactionで検査する。
全scaled/root/products上界1896419、RAW8MB、45秒soft/60秒hard・2GB。
prefixはQ1FHP001一ファイルで、全input/root/accumulatorと全cursorを保存する。
部分prefixからvaluationを主張せず、自動継続しない。出力I/O失敗でfull prefixを重複保存しない。

この新しい単独目的はF valuationだけである。T0、他因子、他tail、GCD、点、
strict cover enlargementや実不良集合の除去を含まない。
