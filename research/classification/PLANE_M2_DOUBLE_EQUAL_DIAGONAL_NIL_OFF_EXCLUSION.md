# equal-diagonal・nilpotent off-field：全higher actualDの因子を強制する

2026-11-08。主義務G/E。新しいpureproofで、source matrix、profile、GB又はjetの追加実行なし。
実Dが既約被約quintic、rank二Poisson plane、normal位数二、constant double rootのh=f型を扱う。
本稿はtangent conic+lineとconcurrent三線の全h=f sourceを全空にする。
全sourceの一部の点又は有限sampleへ限定しない。

## 1. 全source仮定とkの正しいunit

base三座標をb=(x,y,z)、normalを(u,v)とする。
全weight零tensorは

\[
 P_0^{ij}=(P_f)^{ij},\quad P_0^{i,u}=uX_i,\quad
 P_0^{i,v}=uZ_i+vX_i,\quad P_0^{u,v}=k u^2.
\]

X,Zはtracefree linear fields、X(f)=Z(f)=0。
ここでZのmatrixはnilpotent、fはreducible cubicとする。
tangent h=fではX=-2N、Z=cN、N=2y partial_x+z partial_y。
concurrent h=fではX=(ax+by)partial_z、Z=(cx+dy)partial_zである。
両caseでdet(E,X,Z)=0、Z(f)=0なのでactualG2=3f k u^2。
normal位数二により**k nonzeroが結論として従う**。
k=0はnormal位数三以上の別義務で、normalizerで無断に反転したものではない。

全P1はraw36（又は忠実Euler商30）を許し、全P2九を許す。
P1のmixed coefficientはnormal quadraticのconstant倍、P2にはnormal columnsが無い。
従ってfullPのnormal bracketは引き続きk u^2である。
bad u-row coefficient、全off-block及び他の未知higherを零にしない。

## 2. canonical Hamiltonian ratioがbase affine fieldを与える

R_i=sum_j partial_j P^{ij}をrow divergenceとする。
trace X=trace Z=0と全P1/P2の支持から

\[
 R_u=0,\qquad R_v=-2k u.
\]

canonical divergence-zero strict liftPi=P+(E wedge R)/5はactualFに対してPi sharp dF=0。
d(v/u)のEuler contractionは零なので

\[
 M=\Pi^\sharp d(v/u)
   =P^\sharp d(v/u)-\frac{2k}{5}E,\qquad M(F)=0.\tag{1}
\]

このvertical Euler correctionを落とさない。
直接fullbracketsを代入すると

\[
 M(u)=\frac{3k}{5}u,\quad M(v)=\frac{3k}{5}v,
\]

\[
 M(b)=(Z-\tfrac{2k}{5}I)b+B(u,v).\tag{2}
\]

X=Yなので(v/u)(Y-X)項は正確に消える。
Bはbaseに依存しないrational homogeneous degree一の三vectorで、全P1 mixed sourceを含む。
baseP1及び全P2はd(v/u) contractionへ寄与しない。
normal denominator uはdense open上の有理恒等式を読むためだけで、source parameterは反転しない。

## 3. 全rational translationを一度だけ消す

Zはnilpotent、k nonzeroなのでkI-Zは常にinvertibleで、逆は有限nilpotent sumである。
T(u,v)=-(kI-Z)^-1 B(u,v)、b'=b+T(u,v)と置く。
Tはnormal homogeneous degree一であるからM(T)=3kT/5。
式(2)から全translationが消え、

\[
 M(b')=(Z-\tfrac{2k}{5}I)b'.\tag{3}
\]

逆はb=b'-T。これはdense fieldC(u,v)上のaffine invertible changeであり、
global PGL変換と誤認しない。最終因数分解への含意はGauss lemmaで戻す。

## 4. normal gradingの一意性から全Fを復元する

Fをb'で書き、base homogeneous degree dごとに分ける。
全体degree五、Tのnormal degree一から、対応normal degreeはm=5-d、0<=m<=5。
coefficientsはC(v/u)u^mに入る。uとvの個別exponentsが負でもnormal total degreeはmのまま。
式(3)とnormalEulerから各summandに対するMは

\[
 Z+k(m-2)I.
\]

Zは有限degreeのbase polynomial space上nilpotentなのでm not equal twoではinvertible。
ゆえにFの全summandsはm=2以外零。
rational translationはnormal degree一であり、元のnormal order二leadingは変えない。
従って正確に

\[
 F=u^2 f(b')=u^2 f(b+T(u,v))\tag{4}
\]

となる（全体nonzero scaleを除く）。
Z(f)=0によりm=2 summandのinvarianceも整合する。
元P0の未知higherから新しいdegree五項を自由に足す余地はない。

## 5. fraction fieldでの因子はpolynomialFの既約性に反する

K=C(u,v)とし、FをK[x,y,z]で見る。
fはreducibleであり、base-affine invertible translationb to b+Tは各nonconstant factorをnonconstantへ移す。
式(4)はK上で二つ以上のpositive base-degree factorsに分解する。
actualFはC[u,v,x,y,z]上polynomialなので、UFDのGauss lemmaからC上で既約とはなれない。
denominatorのcancellationが起きても、fraction fieldで二positive-degree factorsを持つ結論は変わらない。

従って**tangent conic+line及びconcurrent三線の全constant-double h=f sourceは、
実D既約被約の条件下で空**。
この結論はbadCの方向又はnonzero normalizerを追加に仮定しない全raw36/9 sourceの命題である。

tangent h=f,Z0の旧像22上界は正しいままだが、本稿が受理されれば全空へ強化される。
他のh=L^2M/L^3、moving root、nonreducedR、m3の補集合は別に保持する。
本pureproofの独立読解前にglobal claimを更新しない。
