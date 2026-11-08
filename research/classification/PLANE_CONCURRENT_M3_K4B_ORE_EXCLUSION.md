# K4bの全Ore source：maximal poleと被約cubicの矛盾

2026-11-08、PRIMARY G、producer純証明・独立review前。
新CAS、source probe、matrix、GB又はjetは無い。
全parameter/magnetic/W sourceを保持し、irreducible Dかつreduced concurrent cubicを持つm三源を排除する。

## 1. 全Qと全magnetic正常形

whole K4bは

```
U=(u²,αuv+βv²), V=(0,u²),
A1=0, B1=(2−α)u−2βv,
Qab=bB1+C2, Qan=U,Qbn=V,Quv=0.
```

K4=PfQ=−u⁴。
fullP=Q+∂c∧W、whole moduleは先行lemma通りY=HZ+pi#dH。
irreducible DにはHhomのu-poleが正確に四でなければならない。
全Hはchartv1でΣ_(k=−4)^1 h_k(a,b)u^k、
deg h_-4..h_-1≤二、deg h0≤一、h1 constant。

先行全magnetic正常形により

* β非零：C2零。
* β零,α≠1：C2=εv²、ε0/1。
* β零,α1：C2=0,uv,v²。

を全て保持する。これらの正常化はactual shifts/scalesであり、原C2 sourceの補集合を落とさない。

## 2. β非零はmaximal poleを持てない

chartv1で

```
pi_au=(1−α)u²−βu, pi_bu=−u³,
pi_ab=a u²+2(1−α)bu−3βb,
Z=(αu+β,u²,0).
```

H=h(a,b)u^-4+lower polesとする。
Yuのu^-3係数からβ h_a=0、従ってh=h(b)。
Yaのu^-4係数はβ(5h−3b h_b)。
degree≤二のnonzero polynomialに5h−3b h_b=0を満たすものは無い。
ゆえにh零、Dはu因子を持つ。
β非零の全α/Wを同時に排除した。

## 3. β零、C2零/uvも全actualu因子

d=1−αとし、C2=δuvと書く（δ非零はα1のspecial sourceだけ）。

```
pi_au=d u², pi_bu=−u³,
pi_ab=a u²+2dbu+δu,
Z=(αu,u²,0).
```

d非零ならYuのu^-2からh_a0。
Yuのu^-1から(h_-3)_a=h_b/d。
Ybのu^-2から

`5h−(2b+δ/d)h_b=0`。

finite polynomialのhighest b-degree≤二に対する係数5−2nは非零だからh0。
d零ならYuのu^-1がh_b0を与え、Yaのu^-3係数はh+δh_b=h。
やはりh0。
全α、zero/nonzeroδのactual factor complementを覆った。

## 4. 唯一残るC2=v²：全negative係数の有限再帰

β0,C2=v²、全αを保持する。
chart coefficientsは

```
pi_au=d u², pi_bu=−u³,
pi_ab=1+2dbu+a u²,
Z=(αu,u²,0).
```

全module componentsは

```
Ya=αuH+(1+2dbu+a u²)H_b+d u²H_u,
Yb=u²H−(1+2dbu+a u²)H_a−u³H_u,
Yu=−d u²H_a+u³H_b.                                   (1)
```

Ya/Ybのu^-4から(h_-4)_a=(h_-4)_b0。
従ってh_-4=A constant、A非零。
Yaのu^-3とYbのu^-3から

`h_-3=(4−5α)A b+B`。

次のu^-2 coefficientsは

```
h_-2=5A a+U A b²+V B b+C,
U=−(2α−1)(4−5α)/2, V=3−4α.
```

Yaのu^-1を積分すると

```
h_-1=(6−10α)Aab−(2−α)U A b³/3
      −αV B b²/2−(3α−2)C b+l(a).
```

degree h_-1≤二を戻すと、A非零より

`(2−α)(2α−1)(4−5α)=0`。

したがって全sourceはα=2,1/2,4/5へ強制される。
Ybのu^-1からl'(a)=4B、すなわちl=4Ba+D。
Yuの全negative coefficientsもこのsourceで零となる。
全homogeneous solutionsとlower coefficientsを残しており、generic αだけのkernelへ限定していない。
h0/h1は全linear homogeneous Hhomの四係数である。

## 5. 全remainingHのbase quadratic部分はnonreduced fを強制

T=6−10α、Z=−αV/2=α(4α−3)/2と書く。
上の全Hでquadratic(a,b)部分は正確に

`H_quad=A U b²/u²+(A T ab+B Z b²)/u`。                 (2)

他のA/B/C/D/linear項はbase degree≤一なのでfへ寄与しない。
(1)からY_a/Y_bのbase degree二部分を取り、
W0 base-gradient gaugeへ戻す。
一般にWv0=−(∂aYa+∂bYb)/3だから、全coefficient計算は

```
Wa=K ab+L b²,
Wb=−(K/2)b²,
K=A*((6−4α)T−2U)/3,
L=B*Z*(3−2α).
```

u依存のquadratic(a,b)項はWv0のEuler補正で全て相殺する。
従って元anchorは

`f=(K/2)ab²+(L/3)b³ = b²*((K/2)a+(L/3)b)`。

全α=2,1/2,4/5、全A/B/C/D/linear係数でnonreduced又は零。
stated reduced concurrent cubic hypothesisに反する。
chosen coefficientを零にしたfactorではなく、全maximal-pole sourceに共通する結論である。

## 6. 結論

K4bの全α/β/C2/W補集合は、actualu factor又はnonreduced fで排除される。
よってreduced concurrent R・normalorder三・irreducible quinticの**全K4b源は空**。
先行K211/K31因子排除とK22/K4aの一actual軌道上界を合わせると、
strict Oreの全nonzero K4枝はG像≤24。
残るwhole frontierはK4零の非constant-factor源。
K0、fullE又はambient component数にはまだ未解決がある。
