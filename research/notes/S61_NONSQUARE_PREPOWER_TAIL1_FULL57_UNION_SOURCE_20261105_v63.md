# 二つのtail十分開と次の既知単元valuation source

2026-11-05。独立SOURCE149／CONTROL、tail1独立QQ axis peer、旧DIRECT f946とGLOBAL664を結ぶ。
root peerは元literal14188から149-nodeと四κを独立再構成し、四列の全axisとT1の全axisを照合した。
特にa³e⁰d⁰係数は

```
-211840324830202385382405189416855499140847740295053312
/4503265050648227110358028771123986625 !=0.
```

全7595係数の独立同一性とは言わないが、実式T1の非零性はこれで証明済みである。

## 正規化を含む十五角の恒等式

source順は`[0,1,3,6,2,11,12,9,7,13,10,4,5,8,14]`。
旧DIRECTと同じG18/a、G39/a、各c_i、aG129/c12、G120−G99、行一＋二行二を全列に施す。
source GL15のdetは`-1/(a*c6*c8*c9*c10*c11*c12)`で全FULLU単元。
target順は

```
[22,11,13,25,15,19,14,18,28,30,29,27,31,1,17].
```

H27はL逆、H31は(L gamma)逆、H3三列はL逆、今回のtail1はB倍。
従ってtarget倍率の積は`B/(L^5 gamma)`、`L=aFC²`、`gamma=-4FBC²/625`。
これをD15_1と定義する。新しいκ3やAの反転はない。

旧自由16変数の恒等式`det(A V-c u)=A² detA4`は追加列が任意である。
その列をB倍tail1に替えれば`detA4_1=-T1`となる。
三pair行のanchor消去、正符号のblock Laplace、anchor-first permutationの正符号を保つので

```
A_anchor*D15_1=-T1*D12.
```

元の全sourceでのdomain恒等式であり、Aやκ3の非零chartでしか成り立たない式ではない。
GLOBAL full57橋渡しにより、十分開は正確に

```
FULLU D(D12) intersect (D(T) union D(T1)).
```

旧Omega因子、P、κ3、Aを余分なlocalizerへ戻さない。
T=0で新しい点を覆う、unit ideal、実bad Fitting idealとの等号、全FULLU57は主張しない。

## 保存済み支持だけから分かる有限比較

既知a単元だけを除きU1=T1/a³、U2=T/a⁴とする。保存済み係数の支持は、
U1がa,e,d次数(20,41,18)・全次数58、U2が(20,40,18)・全次数57。
U2がU1を通常のQQ[a,e,d]で割り切るなら商はalpha*e+betaに限られる。
しかしU1の78支持monomialは、U2の支持とそのe倍の支持の合併に入らない。
最初のものは(a,e,d)=(11,0,0)である。新しい積・GCD・factorizationをせず、
**通常の多項式環での割切りだけ**を排除する。
既知FULLU単元の因子や共通の重複成分は未処理なので、局所化後の互いに素や
radical ideal、開の厳密拡大をここから推論しない。

## 次はCという既知単元のvaluationだけを小さく計算できる

C=100(a−a0)、a0=(-18e²+9e+2)/100なので、U(a0+z,e,d)の最初の非零Hasse係数の次数が
厳密にv_C(U)である。両Uのa次数は20。固定整数分母D*100²⁰の下で、各k係数は

```
sum c_ijl * binom(i,k) * 100^(20-i+k)
    *(-18e²+9e+2)^(i-k) *e^j*d^l.
```

一回につき二変数accumulatorだけを持ち、k=0..20を最初の非零まで順に調べる。
三変数商を生成せず、polynomial GCDやfactor searchも行わない。
root numeratorのl1=29<100により、全整数中間はsource input l1×2²⁰×100²⁰以下。
preflightは518 bits以下、owned34188以下、積6528123対以下、RAW32MiB以下。
これは**未実行の次budget提案**である。C valuationの値はまだない。
目的は既知単元による膨張を除くための厳密な高さ情報を得ることで、単独で不良集合を消さない。
