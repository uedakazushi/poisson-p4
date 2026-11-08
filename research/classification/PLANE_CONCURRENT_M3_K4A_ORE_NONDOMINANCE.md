# K4a全Ore源：maximal poleの有限再帰と一actual軌道

2026-11-08、PRIMARY G、producer純証明・独立review前。
新CAS、source probe、module matrix、GB又はjetは無い。
全K4a sourceをnormalorder三/reduced concurrent cubicの下で尽くし、G像上界24を得る。
K4b/K0はまだ別残余である。

## 1. 全K4aと全magnetic補集合

全strict Qは

```
U=(u²,αu²+βuv),
V=(uv,u²+αuv+βv²),
A1=αu+(2β−3)v, B1=(2−β)u,
Qab=aA1+bB1+C2, Qan=U,Qbn=V,Quv=0.
```

全α,βを保ちK4=PfQ=−u⁴。
full P=Q+∂c∧W、Wは全quadratic sourceで、normalorder三とf reducedを元に戻す。
先行whole module lemmaから、chartv1では

`Y=H Z+pi#dH`, D=−u⁴(c+Hhom)。                       (1)

Hhom degree一、u⁴Hhomはpolynomial、normalorder三からbasedegree≤二。
irreducible DにはHhomのu-poleが正確にorder四でなければならない。
lower orderならDはuで割れ、対象外。

β1では全C2 normalformsを戻す。
chartu0のpi_ab=C2(0,1)+O(u)、Z_b=1+O(u)、pi_bu=O(u²)。
H=h(a,b)u^-4+lower polesとすれば、Ybのu^-4係数は

`h−C2(0,1)h_a`。

これは有限polynomials上可逆な1−constant∂aなのでh零。
maximal poleを持てず、**β1の全magnetic sourceはactualu因子で空**。

β≠1ではshift mapはrank三、C2全体をactual positive a/b-shiftsで零へ戻せる。
normal v'=v+k u、base a'=a,b'=b+k aは

`α→α+(1−β)k`, β不変

を実現し、K4a coefficient一を保つ。
従ってk=−α/(1−β)でα0、C2零へactualに正常化出来る。
β1へこのmapを適用せず、上のfactor補集合を保持した。

## 2. 全Laurent sourceとliteral recurrence

d=1−β≠0と置く。chartv1で

```
pi_au=d u², pi_bu=d u−u³,
pi_ab=a(u²−3d)+2d b u,
Z=(βu,u²+β,0).
```

global Hhomのdegreeとnormalorderから

`H=Σ_(k=−4)^1 h_k(a,b)u^k`、

deg h_-4,h_-3,h_-2,h_-1≤二、deg h0≤一、h1 constant。
これが元quinticの全二十八monomialであり、infinite seriesは使わない。
Yの全negative u coefficientsを零にする。

Yuのu^-3,-2,-1から

```
(h_-4)_b=0,
(h_-3)_b=−(h_-4)_a,
(h_-2)_b=−(h_-3)_a.                                  (2)
```

Ya/Ybは

```
Ya=βuH+[a(u²−3d)+2dbu]H_b+d u²H_u,
Yb=(u²+β)H−[a(u²−3d)+2dbu]H_a+(du−u³)H_u.             (3)
```

Ybのu^-4（及び(2)を戻したYaのu^-3）は

`(β−4d)h_-4+3d a(h_-4)_a=0`。

h_-4はaだけのdegree≤二polynomial、nonzeroなので正確に

`h_-4=A a^n`, n=0,1,2、
`β=(4−3n)/(5−3n) ∈{4/5,1/2,2}`, A≠0。               (4)

全βを有限個の任意samplesへ置いたのではなく、irreducibilityのmaximal-poleと元係数から(4)が強制される。

L_j=β+jd+3d a∂aと書く。
(2)(3)の次のcoefficientsは

```
h_-3=−n A a^(n−1)b,                                  (5)
h_-2=n(n−1)A a^(n−2)b²/2−(5−n)A a^n/(2d).            (6)
```

導出はL_-3 h_-3=2d b(h_-4)_a、
L_-2 h_-2=2d b(h_-3)_a−5h_-4+a(h_-4)_aである。
各polynomial homogeneous kernelの必要degreeはn−1/3,n−2/3なので存在しない。
この段でzero/nonzeroのkernel branchを落としていない。

Yaのu^-1から

`(h_-1)_b=−(3n²−15n+10) A a^(n−1)/(6d)`。             (7)

Ybのu^-1が残るhomogeneous coefficientを決める。
これで全negative coefficientsを尽くし、h0/h1はlinear homogeneous Hhomに対応する全四係数として残る。

## 3. n零/一は全m三対象から外れる

n0、β4/5では(7)の右辺がnonzero constant倍のa^-1。
polynomial h_-1は存在せず、このmaximal-pole sourceは空。

n1、β1/2,d1/2では全solutionは

```
h_-4=Aa, h_-3=−Ab, h_-2=−4Aa,
h_-1=(2/3)Ab+B,
h0/h1=全linear homogeneous項。
```

quadratic(a,b) coefficientが無いのでu⁴Hhomはnormalorder四以上、G3零。
今回のnormalorder三に反する。
全B及び全linear Hhomを保持した結論である。

## 4. n二：全module HはAΨ+B av/u+linear

n2、β2,d−1では

```
h_-4=Aa²,
h_-3=−2Aab,
h_-2=A b²+(3/2)A a²,
h_-1=−(4/3)Aab+B a.
```

Ybのu^-1方程式はL_-1 h_-1−2d b(h_-2)_a+4h_-3−a(h_-3)_a0であり、
homogeneous kernelはB aの全一方向。
全global Hhomは

```
Hhom=A Ψ+B*a v/u+ell(a,b,u,v),
Ψ=a²v³/u⁴−2abv²/u³+(b²+(3/2)a²)v/u²−(4/3)ab/u.       (8)
```

全pole条件を元のpolynomial sourceへ戻す必要方向だけでも以下の上界に十分。
module lemmaの一意性から、同じHhomを持つambient tensorはEulerを除いて一意である。
未計算extensionを追加して数えない。

## 5. Canonical anchor、実正常化、全G上界

A部分の全polynomial W0は

```
Wa=(2/9)Aab,
Wb=(3/2)Aa²−(1/9)Ab²,
Wu=−(1/3)Aa v+(2/9)Ab u,
Wv=−(4/9)Ab v+(4/3)Aa u.
```

従ってbase gradient gaugeで

`f=A*(ab²/9−a³/2)`。

reduced fは正確にA≠0であり、A0は先行normalorder/reduced仮定から除かれる。

全ellはactual c'=c+ellで除ける。
β2のQを保つpositive shearはa'=a+t u,b'=b+t v。
Ψについて直接微分すると

`(u∂a+v∂b)Ψ=(5/3)av/u−(4/3)b`。

二次Taylor残りはlinear homogeneous formなので、再びc-shiftへ吸収される。
従ってB→B−5At/3、A非零だから全Bを零へ戻せる。
最後にc-scaleでA一へ直す。全inverseはnonzero scalar又はminus shift。

よって全対象は一actual projective PGL5 orbitへ入る。
canonical Qはα0,β2,C2零、canonical Wは上のA一。
actualDはoverall scalarを除き

`D=u⁴c+a²v³−2abuv²+(b²+(3/2)a²)u²v−(4/3)ab u³`。

c一次、coefficientu⁴と残項gcd一なので既約。
canonical fは三distinct concurrent lines。
元rational module identityからprojective Jacobiはdense openで零、従ってglobal polynomialとして零。
PGL5 orbit及びclosureのdimensionは高々24。
これで**K4a whole G非優越**を得る。
exact orbitdimension、ambient componenthood、既component E所属は未主張。
K4b/K0及びconcurrent G全体はまだ残余である。
