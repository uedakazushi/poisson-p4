# 八HF wall：六全源を排除し、二つの実Poisson A10族へ還元する

2026-11-08。PRIMARY `E_reduced_nonnormal`。producer fullsource算術と純幾何の論証。
新成分性・二つの残族の親・全Eは未証明。旧34/51計算を再実行していない。

## 1. 全sourceと正しい順序

R=2L+M、m二normaldouble、constant primitive rootの全first-symbol正規化は
`PLANE_DOUBLE_LINE_M2_CONSTANT_DOUBLE_ENTRANCE.md` と独立全34reviewが与える。
残る八fixedHF wallは、(k,p)=(0,5/8),(−1/2,2/5),(−2,−1/2),(4,7/4)、各ε=0/1である。

旧reportに保存された**全150 coefficientでの三十basisと九basis**からP1/P2を再構成する。
HFのmixed source順はB,C,D,E（18から各三）で、以前のm三D,B,E,C順を流用しない。
completeK1は外側にkernel columnsを保存するので、転置してfullsourceへ掛ける。
K1次元は8,8,8,8,10,9,9,9、P2 homogeneous核は0,0,2,1,0,0,0,0。
全inhomogeneous P2九を同じ既存rank7/8/9minorで解き、残compatibilityと全homogeneous核係数を戻す。

全P0のu-rowはu倍であり、K1のbad u-rowは全八でC_c=mだけである。
従ってm零なら、**全P2九**にnormal columnsがないためactualDはu因子を持つ。
これはhomogeneousP2だけからの推論でない。

## 2. 六全源の最弱排除

ε一の四源は、元全350への上記inhomogeneous式代入後に次のexact rowを持つ。
source coordinatesを各reportのp_iとし、Jacobiはnative cyclic conventionである。

| k | genuine bad m | 元row | 値 |
|---:|---|---:|---:|
| 0 | p1 | 53 | 256 m²/3 |
| −1/2 | p1 | 53 | 25 m² |
| −2 | p3 | 124 | 16 m²/3 |
| 4 | p2 | 53 | 64 m² |

従ってC上のreduced supportではm零、actualDにu因子がある。
全homogeneous P2 parameterも式に残したうえで、このrowはそれらに依存しない。

ε零のk零源は、全sourceで正確に

`D=5*(-2b+p3*u+p6*v)*(3a*u+8p1*v²+3p2*u²+3p5*u*v)²/48`。

ε零のk四源も全Dに `-6b+p4*u+p7*v` 因子を持つ。
両linear formはb係数が非零定数なので、sourceのzero strataでもproperである。
この六HF sourceは被約既約quinticを持たない。

## 3. 二残源は全parameterでPoissonなA10族

残るε零は

```
A: k=−1/2,p=2/5: K1八＋homogeneousP2二,
B: k=−2,p=−1/2: K1十＋homogeneousP2零.
```

いずれもinhomogeneous P2を戻した後、**元全350Jacobiが全parameter polynomialとして恒等零**である。
従って単なるnecessary dimension sourceではなく、明示した全A10がactualPoisson familyを与える。
各familyのPi係数はparameterのpolynomialで、P0のplane restrictionは常に同じ非零tensorなので、
射影Poisson spaceへのmapは全A10で定義される。
Bのc依存Pab termはそのまま残す。全族がc-Oreと仮定しない。

旧独立34reviewのwhole image上界は両HF sourceについて23。
これは固定kernelの点だけでなく、全normalizer inverse、plane、PGL作用を戻した全source像の上界である。
新A10がそのwhole実sourceの中にあることから、各PGL5飽和閉包Y_A,Y_Bは
**既約でdimension≤23**となる。

## 4. 最弱のdense-orbit bridge：全transitivityは不要

一般にirreducible parameter family TからのactualPoisson像をPGL5で飽和した閉包Yが既約、dimY≤23とする。
同じfamilyの一つのactual点Pの射影PGL軌道が23次元なら、
その閉包はYに含まれる23次元既約閉集合なので **Y=closure(PGL5·P)**。
従ってTの全点がこの一つの軌道閉包に入る。全Tで同じorbitであることを証明する必要はない。

さらにその**同じ点**のfullambientDがrank102、Euler24/実M51rank48なら、
raw局所次元48が接次元48と一致し、Yはsmooth unique23の全ambient componentである。
D101ならこの全familyのorbitclosure23がより大きい親の境界であり得るため、親のactual認識を別に行う。
D≤100なら局所profile列へ進まず、actualparent幾何に戻る。
このbridgeは一般fiberのflatness、全T上のorbit次元外挿、sourceclosure上の任意fiberを仮定しない。

## 5. familyから構造的に固定できる二anchor

未実行の後続成分性判定のため、各familyのirreducibleDを持つ一点を固定した。
点列又はrank/profile探索ではなく、bad mとpolar coefficientを一にしてprimitive-linear UFD証明を与える。

Aはp1=q1=1、他のhigher coefficients零。primitive actualDは

`F_A=25c*u⁴−4b*(3a*u+5v²)²−20a*u³v−50u²v³`。

c係数25u⁴と残部は互いに素（bv⁴係数が非零）なのでC上絶対既約。
同じtensorを整数化するscaleは20。

Bはp5=p2=1、他higher coefficients零。primitive actualDは

`F_B=−b*(3a*u+2v²)²−6a*c*u²v+2a*v⁴−12c²u³+4c*u*v³`。

b-linearであり、係数quadratic-squareのprime quadraticはu非可除、
残部のc²係数は非零定数*u³なので、そのquadraticが残部を割らない。
従ってprimitive-linear UFDからC上絶対既約。同じtensorの整数scaleは2。

二つのfullPi、実Pfaffian、primitive比、全parameter選択はanchor receiptに保存した。
全族のJacobi恒等零からこれらもPoissonである。
**まだこの二点のambientD/orbit matrixは計算していない。**
一つずつ同じ点だけのfullD/実Mが§4の最弱後続判定になり得る。
新component数又は既知親を本稿では推測しない。

## 実行と保存

一回の許可された全八source batch、gate PASS後180秒/3GiB。
9.433秒/RSS133640KiB、全算術PASS。全P1/FULL9P2+homq、全350、fullD/Pf、badm、
実Euler gaugeの成功／不成功をreportに保存した。
shared6022は編集せず、exactinlinecallerをjob/archiveへ固定した。
原source/全壁/全逆の独立受理と、この新source polynomial恒等式／新排除の独立確認は区別する。
