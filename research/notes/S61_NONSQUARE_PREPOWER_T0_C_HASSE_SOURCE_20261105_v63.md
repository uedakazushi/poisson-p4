# 旧tail T0/root239のC-Hasse SOURCE予算

2026-11-05。T0は既存T=Bκ*tail2、KAPPA_EXTRA2のroot239である。
7208係数の別の点や新しいtailを評価するものではない。
旧全QQ source peerとliteral hashを保持し、a⁴というFULLU単元を除いたU0を入力とする。
今回、C冪・Hasse層・C valuationは未計算で、workerもまだ構築していない。

**offsetの規約を固定する。** root239.polyの最初の4はa-offset、次が7208項数。
bodyのa指数は既に0..20であり、そのままU0=T0/a⁴の指数である。
original .jsonのa指数は4..24なので、そちらを使う場合だけ4を引く。
二重にshiftしない。元のliteral順序も変えない。

保存済み共通D197へ定数だけを直すと、integer coefficient最大320 bits、l1は326 bits。
`d|D, scale=D/d, scale*d=D, I=n*scale`の小さな検査だけを使い、巨大なQQ cross-productは作らない。
新しいinteger coefficientファイルはまだ作っていない。

C=100a+18e²−9e−2、a0=(2+9e−18e²)/100。
a次数20、e次数40、d次数18なので、k0..20を最初の非零までstreamする。
root l1=29<100、input_l1<2³²⁶から、全integer partialは481 bits以下、hard518を維持できる。
root cache441、accumulator e≤80,d≤18の1539、input7208とcoef/key双方、metadata16000でowned34376。
全conversion/scaled/root/products上界4117089、RAW6.582MB、hard8MB、45秒soft/60秒hard・2GBの提案である。

前準備もplain sourceを一行ずつ読み、full .jsonとinteger payloadを同時に抱えない。
fresh source/code独立監査と親の一回GO以前にCの値を出さない。
T0/T1を小さくするための入力整備であり、GCD、他因子、不良集合や補集合の除去は含まない。
