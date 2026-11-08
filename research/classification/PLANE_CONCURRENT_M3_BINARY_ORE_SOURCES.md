# Strict binary Ore sourceの全actual正常形

2026-11-08、PRIMARY G、producer純証明・独立review前。
新source probe、minor、GB、jet又はCASは無い。
先行whole Ore bridgeからの全strict Qを、実block GL2×GL2とtensor scaleで有限個の低parameter族へ移す。
全central coefficient三つ及びfull W compatibilityを保持する。
本稿だけで全concurrent G又はEが閉じるとは主張しない。

## 1. 全sourceと許される実作用

y=(a,b,u,v)、normal planeのEuler E_n=u∂u+v∂v。
全sourceは

```
Q=(a A1(u,v)+b B1(u,v)+C2(u,v))∂a∧∂b
  +∂a∧U2(u,v)+∂b∧V2(u,v),
[U2,V2]=A1 U2+B1 V2,                                  (1)
C2=c0*u²+c1*u*v+c2*v².
```

全ΠはQ+∂c∧W(y)、全quadratic Wに元のnormalweight/first-symbol条件を戻し、
exact moduleはE_y∧L_WQ=0である。
Qを正常化するため(a,b)へ全constant GL2、(u,v)へ全constant GL2を使える。
これらはcを固定するactual PGL5であり、元RのmarksやanchorはWへ戻す。
Qのtensor common scaleもfull Πのscaleとして戻す。
base a/bを変えた後に元のbinary cubic fが同じliteral式のままとは仮定しない。
full Wを含むcoverであり、変換後の元first-symbol/leading conditionsを必要なら明示的に課す。

K=det(U2,V2)とする。Kのnonzero scalarはnormalizeしてよいが、K0を別sourceとして保持する。
K≠0では(1)のA1/B1は一意。
全normalizerの逆はconstant invertible block matrixとnonzero scalarの逆で、birational変換を使わない。

## 2. K非零・三roots（211）の全一parameter族

先行whole bridgeによりKは211/31/22/4だけ。
三distinct rootsをu,v,u−vとしg=u v(u−v)。
全quadratic logarithmic fieldはl1 E_n+t0 J_g、J_g=(g_v,−g_u)。
pairをbaseGL2で動かし、J_g係数を第一row零/第二row一に出来る。
K非零より第一rowはnonzero linear l1 E_nで、detは−3l1*g。
第四rootを持てないためl1はdouble line。これをuとする。
row scalesでU=uE_n、V=J_g+l2E_nとし、V→V+sUでl2のu成分を除く。
従って全sourceは

```
U=(u²,u v),
V=(u²+(t−2)u v, −2u v+(t+1)v²),
A1=−u+(2−t)v, B1=u,       t∈C.                       (K211)
```

K=−3u²v(u−v)。全tで(1)は直接のLie bracket恒等式。
marked double lineの全三choices及びnormal/base framesは上のactual作用で戻る。

## 3. K非零・二rootsのmatrix-pencil全作用

二rootsをu,vへ取り、全logarithmic fieldsを

`U=(u l1,v m1), V=(u l2,v m2)`

と書く。L=[[l1,m1],[l2,m2]]=uA+vBは全linear pencil。
base GL2はLの**左**に作用する。
Euler directionsを任意に混ぜる右GL2をactual変換だと仮定しない。
normal diagonalではu,vのargument scalesだけを使う。

### 3.1 detL∝uv（K22）

A,Bは各rank一、image lines及びrow covectorsは独立。
左GL2でimage linesを座標軸へ取ると

```
U=(α1 u², α2 u v),
V=(β1 u v, β2 v²),
A1=−β1 v, B1=α2 u.                                   (K22)
```

α=(α1,α2)、β=(β1,β2)はordered P1 pair、det(α,β)≠0。
各row scaleがprojectivizationを実現する。
全四standard affine chartsを使い、zero entryを除かない。
高々二parametersで、K=det(α,β)u²v²。
この左右の区別により、α/βを任意conjugacyで一つのpointへ潰さない。

### 3.2 detL∝u²（K31）

Aはinvertible、左GL2でA=Iへ直せる。
det(uI+vC)=u²からtrC=detC=0。
normal u/v scalesと左common scaleでCをnonzero constant倍出来る。
従ってCは零、又はprojective nilpotent coneの全P1。

`C=[[r,s],[t,−r]], r²+st=0` に対して

```
U=(u²+r u v, s v²),
V=(t u v, u v−r v²),
A1=(r−2t)v, B1=u+(2r−s)v.                             (K31)
```

全P1は(r,s,t)=(z,1,−z²)、z∈C、及び(0,0,1)で覆う。
C零も別に保持する。K=u³v。
これはnilpotent Cをconjugacyで一つにする操作ではなく、**scalarだけ**を除いた全conic。
K∝uv³はnormal swapで同じ全coverへ入る。

## 4. K非零・一root（K4）の全二族

rootをu0へ取り、U=(u l1,Q2)、V=(u l2,R2)。
detからl1R2−l2Q2=k u³、k≠0。
span(l1,l2)のrank零は不可能。

rank二なら左GL2で(l1,l2)=(u,v)。
uR2−vQ2=k u³からQ2=uL1,R2=k u²+vL1。
normal/base scalesでk=1へ直し、L1=αu+βvとすると

```
U=(u², αu²+βu v),
V=(u v, u²+αu v+βv²),
A1=αu+(2β−3)v, B1=(2−β)u.                            (K4a)
```

全α,β∈Cを保つ。
必要なnonzero square root scaleは複素体上のactual diagonalであり、k0へ適用しない。

rank一なら左GL2でl2零。
l1 R2=k u³からl1∝u,R2∝u²であり、row scalesで
U=(u²,Q2),V=(0,u²)へ直せる。
row shear U→U+sVでQ2のu²係数を除き

```
U=(u², αu v+βv²), V=(0,u²),
A1=0, B1=(2−α)u−2βv.                                 (K4b)
```

全α,βを保持する。
両族K=u⁴で、rank一/二と全zero parameterを同時に覆う。

## 5. K零：primitive共通directionの全次数

K0、(U,V)≠(0,0)とする。
C[u,v]のUFDにより、primitive vector field R_eを共通directionとして

`U=P_(2−e) R_e, V=Q_(2−e) R_e`, e=0,1,2

と書ける。R_eの二componentsのgcdは一。
P,Qがconstant proportionalならeに関係なく次節へ入る。
そうでなければe二は不可能で、e零又は一だけが残る。
factorizationをfunction fieldのまま使わず、primitive polynomial pairへ戻した全次数である。

### 5.1 constant proportional pairは全essential三

U=λR,V=μR、(λ,μ)≠(0,0)なら(1)からλA1+μB1=0。
v0=λ∂a+μ∂bと置くとv0∧Q=0、L_v0 Q=0。
先行独立 `PLANE_DEPENDENT_ORE_SEPARATE_REVIEW.md` §3のconstant-factor補題により、
full module式E_y∧[Q,W]=0を満たす任意WでactualDはcとv0の二方向に依存しない。
既約非零quinticならessential三であり、既受理全八成分/global上界23へ入る。
normal/base frames、planeの選択次元をglobal上界へ再加算しない。
全Rのnonzero/zero componentを含み、chosen constant-factor pointの議論ではない。

U=V0ならQ=F2∂a∧∂b。
F2非零ならactualD=F2 times homogeneous cubicでproper factorを持ち、
F2零ならactualD零。従ってこの全sourceは既約quinticの対象外。

## 6. K零・constant primitive direction：全四pencil/marked-point類

R0をnonzero constant normal field、U=P2R0,V=Q2R0とする。
P2,Q2がindependentの場合、そのpencilはbinary quadraticsの二次元subspace。
dual binary quadraticのrank一ならcommon-root pencil、rank二ならbasepointfree pencilで尽くされる。

common-root pencilのnormalizerはそのpointを保つBorel。
R0のdirectionはそのpoint又はcomplementの二orbits。
basepointfree pencilのnormalizerは二ramification pointsを保つtorusとswap。
R0のdirectionはramification point又はそのcomplementの二orbits。
いずれもmarked R directionを忘れずに分類した。
全normal/base GL2のactual作用から以下の四anchorとなる。

| 型 | U | V | A1 | B1 |
|---|---|---|---|---|
| C_tan | u²∂v | uv∂v | u+t v | −t u |
| C_trans | u²∂u | uv∂u | (1+t)v | −(2+t)u |
| B_axis | u²∂u | v²∂u | 0 | −2u |
| B_other | u²(∂u+∂v) | v²(∂u+∂v) | 2v | −2u |

common-root二類では全t∈Cを保つ。
一般にparticular A1=R0(Q2),B1=−R0(P2)であり、
gcd(P2,Q2)=1のbasepointfree類ではdegree一syzygyが無く一意。
common-root P2=L L1,Q2=L L2ではsyzygyはt(L2,−L1)全て。
これで表の全自由parameterと全zero補集合が導かれる。

## 7. K零・primitive linear direction：全SS二parameter/Jordan一parameter族

R1のmatrixがrank一なら二componentsに共通linear factorがあるのでprimitive性に反する。
従ってR1はinvertible linear normal field。
independent P1,Q1をbaseGL2で(u,v)へ取り直し、U=uR1,V=vR1。
normal GL2とpaired base inverseでR1をJordan正常形へ直す。
eigenvalue common scaleは全Q/coordinate scalesで除ける。

全semisimple（scalar含む）はR1=u∂u+t v∂v、t≠0で

```
U=(u²,t u v), V=(u v,t v²),
A1=(t+s)v, B1=−(1+s)u,       t≠0, s∈C.              (Z_SS)
```

全invertible mixed JordanはR1=u∂u+(u+v)∂vで

```
U=(u²,u²+u v), V=(u v,u v+v²),
A1=u+(1+s)v, B1=−(1+s)u,       s∈C.                  (Z_J)
```

一般式A1=R1(v)+s v,B1=−R1(u)−s uから全syzygyを保持した。
t1のscalar Rも除かない。
t0のrank一fieldはC_transに、nilpotent zero-eigenvalue境界はC_tanに既に入る。
従ってinvertibilityを仮定したchartのcomplementも全K0 sourceに含まれる。

## 8. 全C2と次のwhole module義務

全表でC2の三coefficientsは自由に残す。
actual positive shifts a→a+ell_a(u,v),b→b+ell_b(u,v)はU,V,A1,B1を変えず

`C2→C2+U(ell_b)−V(ell_a)−A1 ell_a−B1 ell_b`。          (6)

各ellは全binary linearで、全inverseはminus shift。
この小さいlinear mapのimageを未計算のままC2を零に固定しない。
future normalizationで使う場合はrank dropの全parameter lociと残magnetic coefficientsを保持する。

以上でstrict Qは、K非零の全五種類（K22は全四charts、K31はP1＋zero）と、
K0の四anchor/SS/Jordan、global essential三又はfactor枝へ全て入る。
各unclosed canonical familyは高々二parametersとC2全三coefficientsを持つ。
whole W moduleのfull sourceは全quadratic W（必要なら元graded sourceのsupersource）で、
E_y∧L_WQ0及び元のactualD3/R条件を戻して判定する。
これらがまだ無い段階でfixed module次元や全G≤24を断言しない。
先行SFの330candidate有限constructionは未実列挙のまま保持し、
本稿のsupersource幾何で将来Gが閉じた場合だけ、その列挙がGには不要だったと結論出来る。
