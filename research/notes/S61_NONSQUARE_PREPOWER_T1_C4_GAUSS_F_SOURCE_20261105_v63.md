# C四乗商の強制scalar正規化と次のF-Hasse SOURCE

2026-11-05。新しい商、factor/GCD、F valuationは計算しない。
保存済み4811係数の整数bitとl1だけを読むSOURCEである。
実scalar正規化やF workerには、C四乗商の独立peerと新しい親budgetを要する。

元の整数表示で `C⁴*Qint=100²⁰*I`、I=(T1/a³)_integerと分かった。
Cの整数係数(100,18,−9,−2)は`−1*(−9)+4*(−2)=1`によりprimitive。
Gaussのcontent補題でC⁴もprimitiveなので、**Qintの全係数は100²⁰で割り切れる**。
これはcontent探索ではなく、固定された単一整数による強制的な整除である。

従って同じQを `K/D1`、`K=Qint/100²⁰`と表せる。
denominatorは人工的な330 bitsから元の197 bitsへ戻る。
まだ新しいK係数ファイルを作っていない。将来のcopy時に各integerの剰余零と逆積だけを確認し、
大きなQQ cross-productは行わない。
raw Qintのcoef428 bits、l1 433 bitsと100²⁰≥2¹³²から、Kのcoef296 bits以下、l1<2³⁰¹が従う。

この正規化後、次の既知単元F=500a+90e²+3e−24に対するa0は
`(24−3e−90e²)/500`。root numeratorのl1は117<500。
a次数16なので、Hasse係数k=0..16を最初の非零までstreamすればよい。
各integer partial sumはinput_l1*2¹⁶*500¹⁶以下で463 bits。
固定scalar正規化の逆積も429 bits以下なので、**hard518 bitsを上げる必要はない**。

input4811、root cache289、accumulator e≤65,d≤18の1254、coefとkey双方、metadata16000により
owned28708。全scaled/root/productは190万以下、RAW8MB、45秒soft/60秒hard・2GBという
新しい単独目的のpreflightである。F冪・Hasse係数・他tailはまだ一切評価していない。
このSOURCEはFULLUの不良集合やT=0の補集合を消すものではない。
