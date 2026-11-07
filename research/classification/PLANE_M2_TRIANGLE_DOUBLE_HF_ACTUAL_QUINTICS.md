# triangle・constant double root h=f：全actualDの二正常形

2026-11-07。主義務G/E。実D被約既約quintic、rank二Poisson plane、zero cubic f=xyz、
normal位数二、constant double root G2=3f u^2を仮定する。
別leading削減のh=f型を扱う。h=L^2M、concurrent、nonreducedR、m3は別残余。
本稿は全actualDの必要正常形を手証明し、新しいlinear profile、GB又はjetは計算しない。

## 1. 全sourceと正しいunit

`PLANE_M2_DOUBLE_ROOT_TRIANGLE_LEADING_REDUCTION.md` の訂正済み全源より

\[
 X=Y=(0,-y,z),\qquad Z=(s x,t y,(-s-t)z),
 \quad P^{u,v}_0=k u^2,\quad k=1+s.
\]

したがってraw kは零でもよい。unitなのはleading amplitude k-s=1であってkそのものではない。
全P1のbad u-rowはc v^2 partial_x、c nonzeroが既約actualDに必要。
他の全30/9係数を保持する。
特にfullPのu-columnは

\[
 P^{x,u}=b_xu^2+a_xuv+cv^2,\quad
 P^{y,u}=-yu+b_yu^2+a_yuv,\quad
 P^{z,u}=zu+b_zu^2+a_zuv,\quad
 P^{v,u}=-k u^2.
\]

係数a_i,b_iはconstantで、normal bracketは全filtrationを通じてk u^2のままである。
既忠実30 Euler商gaugeではb_iを消せるが、ここではraw36の追加u^2項まで許し、
v-columnの追加v^2項も許す。以下のspan(u,v)結論はこれらを含めて不変。
raw fullPをstrict Poissonと仮定しない。

## 2. canonical strict liftの全divergenceを読む

R_i=sum_j partial_j P^{ij}（row divergence）とする。
Rはlinear、div R=0。全P1 mixed columnsのbase dependence、全P2の支持を直接使うと

\[
 R_u=0,\qquad R_v=-2k u,\quad R_x\in\operatorname{span}(u,v),
\]

\[
 R_y=-2y+\operatorname{span}(u,v),\qquad
 R_z=2z+\operatorname{span}(u,v).\tag{1}
\]

P_fのbase divergenceは零、trace X=trace Y=trace Z=0。
P1のbase bracketをbase変数で微分すればu/vのconstant倍だけが残る。
P1 mixedをnormal変数で微分しても同じ。P2のbase derivativeは零なので、式(1)は全sourceで正確。

canonical strict liftはPi=P+(E wedge R)/5。
この符号ではrow divergence(E wedge R)=-5Rとなり、Piはdivergence零。
既受理のcanonical-lift恒等式からactual FはPi sharp dF=0を満たす。
従ってdense u nonzero上のderivationXi=(Pi sharp du)/uはFを殺す。
これはglobalEuler liftをstabilizerと誤認する主張ではなく、canonical strict Casimirの利用である。

R_u=0なのでXiはP-column/u-R/5である。式(1)から

\[
 \Xi_u=0,\quad\Xi_v=-3k u/5,
 \quad\Xi_x=cv^2/u+\operatorname{span}(u,v),
\]

\[
 \Xi_y=-3y/5+\operatorname{span}(u,v),\qquad
 \Xi_z=3z/5+\operatorname{span}(u,v).\tag{2}
\]

actual linear shears y/z to y/z+linear(u,v)で後二式のoff-blockを消せる。
そのSylvester equationsはeigenvalues plus/minus3/5とnilpotent(u,v) blockの間なので常に可逆。
unknown parameterをunitと仮定せず、固定nonzero数だけで解ける。

## 3. semisimpleとlocally nilpotentを分離する

Xiを-5/3倍して記号を改めると

\[
 \widetilde\Xi=S+N,\quad S=y\partial_y-z\partial_z,
 \quad N(u)=0,\quad N(v)=k u,
\]

\[
 N(x)=a v^2/u+bv+d u,\qquad a\ne0.\tag{3}
\]

a nonzeroはbad coefficient c nonzeroそのもの。
NはC[u,u^-1,v,x]上locally nilpotent。
例えばweight(x)=3、weight(v)=1、weight(u)=0のfiltrationをstrictly下げる。
SとNはcommuteする。
Fをy/z weightごとに分けると、nonzero weightでS+Nは有限N-orbit上可逆なのでその成分は零。
従ってFはC[u,v,x,yz]に属す。

k=0ならN=(a v^2/u+bv+d u)partial_xであり、coefficientはa nonzeroのため非零。
N(F)=0はpartial_x F=0を強制する。
しかしS沿いnormal order二の項は3xyz u^2でありxに非自明に依存するので矛盾。
従って**指定の全sourceでk nonzeroが結論として従う**。
k=0を最初から反転していない。

## 4. 全invariant ringとprimitive cubic J

k nonzeroなら

\[
 J_0=u^2x-\frac{a}{3k}v^3-\frac{b}{2k}uv^2-\frac{d}{k}u^2v
\]

はN-invariant。v to v+lambda uでuv^2項を消し、
x to x+mu v+nu uでu^2v及びu^3項を消す。
a/k nonzeroの尺度とcomplex cube rootでv/xをrescaleし、

\[
 J=u^2x-v^3
\]

へ戻せる。全てactual invertible linear transformationsで、Sを保つ。
各段の逆は逆shearと逆nonzero scaleで与えられる。
u以外のrational coordinateをglobal projective equivalenceだとは呼ばない。

u nonzero上はC[u,u^-1,v,x]=C[u,u^-1,J,v]、N=k u partial_vなので
N-invariantsはC[u,u^-1,J]。T=yzを加えると

\[
 C[u,u^{-1},J,T]\cap C[x,y,z,u,v]=C[u,J,T].\tag{4}
\]

実際uの最小negative powerがあれば、そのnumeratorをu零に還元して
P(J,T) to P(-v^3,yz)が零になる必要がある。
J,Tのこの代入はinjectiveなので、最小powerのnonzero numeratorとは両立しない。
これでu-adic cancellationを全parameterで排除し、normalization/base-changeの一般可換性などを使わない。

## 5. 全degree五候補とnormal-order条件

deg(u,J,T)=(1,3,2)なので、全Fは

\[
 F=A u^5+B u^2J+C u^3T+D JT+E uT^2.
\]

uT^2はnormal order一なのでE=0。
全線形shearsはSを保ち、leading F2=unit xyz u^2を変えないためD nonzero。
x to x+(C/D)uと置いてJ'=J+(C/D)u^3にすれば
u^3T項は消え、その逆操作もactuallinearである。
全体nonzero scaleを吸収して

\[
 F=(u^2x-v^3)(yz+B_0u^2)+\kappa u^5.\tag{5}
\]

kappa=0なら二因子に分解するので、既約actualD部分ではkappa nonzero。
B0=0を捨てず、B0 nonzeroとの二whole枝に分ける。

B0=0ではold z=kappa new z、他は恒等としてF=kappa F0。
B0 nonzeroではold u=U、old x=(kappa/B0)X、old v=r V、old y=Y、old z=B0 Z、
r^3=kappa/B0とするとF=kappa F1。
両方とも全逆linear mapが存在し、complex cube rootの有限選択は全点を覆う。
従って**全actualDは次の二projective正常形だけ**へ還元される。

\[
 F_\varepsilon=(u^2x-v^3)(yz+\varepsilon u^2)+u^5,
 \qquad\varepsilon=0,1.\tag{6}
\]

両Fは複素数上で絶対既約である。
xに関する一次式のcoefficientはu^2(yz+epsilon u^2)、constant termは-v^3(yz+epsilon u^2)+u^5。
共通因子uは無く、yz+epsilon u^2のどの因子もconstant termを割らない。
従って係数gcdは一、primitiveな一次式としてGauss lemmaから絶対既約となる。

## 6. 必要最弱の次のfinite source命題

各Fのcanonical strict quadratic sourceは、全150係数に対するdivergence零かつPi sharp dF=0のlinear space LF。
Poisson Jacobi及びactualPf=unit Fを課す前の必要supersourceとして使える。
F1はtorus y to t y,z to t^-1 zを持つためactualD orbit dimension23以下。
F0は更にweights(x,y,z,u,v)=(-2,5,0,1,0)のtorusでweight五なのでorbit22以下。
二torusはprojectively独立である。
従ってdim LF0<=3かつdim LF1<=2なら、全projective源はそれぞれ
22+3-1、23+2-1以下となりG24が閉じる。
fullLFの全Pfaffianが零なら、その型はrank四全空という更に弱い証明でも十分。

これは二wholeactualDの必要sourceを一回有限に解く正当な入口であり、
新しいpoint profile又は旧scalar/m4計算の反復ではない。
現時点では本手証明とLF次元判定は独立未確認、G/E完了は未主張。
