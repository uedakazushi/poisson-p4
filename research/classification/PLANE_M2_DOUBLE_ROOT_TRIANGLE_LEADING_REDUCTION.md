# nonconcurrent triangle・m2 double rootの全leading削減

2026-11-07。主義務G/E。実D既約被約、rank二Poisson plane、zero cubic f=xyz、
S沿い位数二、normal quadratic generic doubleを仮定する。
moving root degree一は `PLANE_M2_MOVING_DOUBLE_ROOT_NORMAL_FORM.md` で全空。
本稿はconstant rootの全P0を四つの低parameter源へ還元する。
全higher帰属又はG次元上界は未完了。計算jobなし。

## 1. 全raw weight零源と全bad-row必要式

G2=gamma h u^2、gamma nonzero、root degree零。
f固定後、hの全component exponentsは111、21、3であり、permutationを戻すと全marked choicesを覆う。
canonical tracefree raw sourceは

\[
 P^{ij}_0=(P_f)^{ij},\quad P^{i,u}_0=uX_i,\quad
 P^{i,v}_0=uZ_i+vY_i,\quad P^{u,v}_0=k_{20}u^2+\kappa uv.
\]

X,Yはtracefree diagonal fieldsでX(f)=Y(f)=0。
Delta=X-Yはconnection relationによって固定される。
actual Euler-Pfaffianの二係数は

\[
 [G_2]_{uv}=3f\kappa-\det(E,X,Y)=0,\qquad
 [G_2]_{u^2}=3f k_{20}-\det(E,X,Z)-Z(f)=\gamma h.\tag{1}
\]

全off-diagonal Zをここで保持する。
fullP1のbad u-row係数v^2 CはC nonzeroが必要。
projective Jacobiのweight一とcurlは

\[
 E\times(\kappa C+DY C)-(2Y-X)\times C+\nabla(C(f))=0,\qquad
 (2\kappa I+4DY-DX)C=0.\tag{2}
\]

高次の他の係数はこの必要式へ寄与せず、C=0ならactualDにu因子が生じる。

## 2. 全closed rational off-diagonal module

omega=(f/h)bはclosed、omega=N2/h、N2=f bはdegree二horizontal one-form。
raw fieldはiE iZ volume=-3h omegaで戻る。

* h=fではomega=bはclosed simple-poleなので全logarithmic二dimensional源。
  Zは全tracefree diagonal fieldである。
* h=x^2yでは、全omegaはd(A y/x+B z/x)+c dlog(x/y)。
  従って全ZはA' y partial_z+B' Delta+C' x partial_z、Delta=diag(-1,2,-1)。
* h=x^3では、全omegaはd(P2/x^2)、P2は全quadratic modulo C x^2。
  従って全Zはspan{x partial_y,z partial_y,x partial_z,y partial_z,diag(0,3,-3)}。

closed-form記述の完全性は初等的に確認出来る。
x=1 chartでh=x^3のomegaはpolynomial closed one-form、係数degree二以下。
horizontal homogeneous numeratorのdegree二条件はpotentialのdegree三部分を消す。
よってpotentialはdegree二以下で、constantを除いた全五係数が上のP2に等しい。
h=x^2yではr=y/x、s=z/xにより高々simple r-pole。
closednessからconstant residueをc dlog rとして引くとpolynomial closed formになる。
同じhorizontal numerator条件がpotentialのdegree二部分を消すので全potentialはAr+Bs。
これでnonlog off-diagonalやzero coefficient境界を捨てていない。

## 3. h=L^3は全P0が空

h=x^3の場合、全Zのx-componentは零、YとXはdiagonal。
従って式(1)のu^2係数はideal(y,z)に属する。
3f k20はもちろんそのidealに属し、det(E,X,Z)とZ(f)も各項にy又はzが残る。
target gamma x^3、gamma nonzeroとは両立しない。
permutationを戻すと全h=L^3 sourceが空であり、higher条件を計算する必要はない。

## 4. h=f：有限diagonal型と全二parameter off-field

Delta=0からX=Y、kappa=0。
curl式は3DY C=0である。
Y=0ならC(f)=0からC=0。
Y nonzeroではC nonzeroはYのzero eigenspaceに属し、trace零よりzero eigenvalueは一つだけ。
座標permutationでC=c partial_x、c nonzeroとする。
式(2)の残係数はY=diag(0,-1,1)を強制する。
全off-field Zは任意のtracefree diagonal二parameterを保つ。
Z=diag(a,b,-a-b)ではdet(E,X,Z)=3a fなので、式(1)はk20-a nonzeroを要求する。
normal u/vの相対scaleはZとk20を同時にscaleするため、k20-a=1へ戻せる。
単にk20 nonzeroとは置かず、k20=0かつa nonzeroの全境界も保つ。
従ってこの全P0は

\[
 X=Y=\operatorname{diag}(0,-1,1),\quad
 Z=\operatorname{diag}(a,b,-a-b),\quad k_{20}=1+a,\quad\kappa=0
\]

の全二parameter源とそのactual linear飽和で覆われる。
normal shear v to v+c uはこの場合Zもk20も変えず、隠れたparameter除去には使わない。

## 5. h=x^2y：全Pf条件と三つのbad-row層

Delta=diag(-1,2,-1)、Y=diag(s1,s2,s3)、sum si=0とする。
式(1)からkappa=s3-s1。
全Z=A y partial_z+B Delta+C0 x partial_zを代入すると

\[
 G_2/u^2=3f(k_{20}+B\kappa)
 -xy(Ay+C_0x)(s_2-s_1+4).
\]

従って

\[
 R=s_2-s_1+4\ne0,\quad A=0,\quad C_0\ne0,\quad k_{20}=-B\kappa.
\]

normal shear v to v-BuでZのB Deltaを消すとk20も零になる。
相対normal scaleでC0を一へ戻し、全P0はZ=x partial_z、k20=0に規格化出来る。

式(2)はCの三座標ごとに異なるbase monomialsを持つため、各nonzero componentは別々に次の二係数を零にする。
di=Deltaのdiagonal entriesとして

\[
 \begin{array}{c|cc}
 C_x\ne0&\kappa+s_1-s_3+d_3+1&\kappa+s_1-s_2+d_2-1\\
 C_y\ne0&\kappa+s_2-s_1+d_1+1&\kappa+s_2-s_3+d_3-1\\
 C_z\ne0&\kappa+s_3-s_2+d_2+1&\kappa+s_3-s_1+d_1-1
 \end{array}
\]

全解は次の三つだけで、互いのY-lociは交わらない。

| bad C direction | 全Y | kappa | nonzero補集合 |
|---|---|---|---|
| partial_x | diag(-2s-1,s+1,s) | 3s+1 | s not equal -2 |
| partial_y | diag(0,2,-2) | -2 | R=6 |
| partial_z | diag(-2,3,-1) | 1 | R=9 |

第一行のR=3(s+2)であり、s=-2ではleading amplitudeが零になるためm三以上の別義務。
この補集合を全空又は他のfamily同値とは扱わない。
Cのmagnitude、他の全P1及び全P2はなお未知sourceとして保持する。

## 6. 次の有限whole-source義務

triangleのmoving rootとconstant h=L^3は全空。
constant h=fの二parameter源とh=L^2Mの一parameter/二fixed源、合計四源が残る。
表示は全normalizationのactual inverse、全部零/非零補集合を保つleading/source reductionである。
これだけで全sourceを既知成分へ含めたり、kernel粗上界をG boundと呼んだりしない。
次のfinite jobが必要なら、この四whole源に対する最弱fullhigher命題と全parameter wallsを先に契約する。
concurrent三線、tangentconic、nonreducedR、m3等は別に保持する。
