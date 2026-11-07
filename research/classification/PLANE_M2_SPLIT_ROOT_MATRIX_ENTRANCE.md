# rank二plane・m2 squarefree normal quadraticの全root-matrix入口

2026-11-07。主義務G、指定一枝ではEも被覆する数学的還元。
S=P2はrestriction rank二、zero cubic f≠0、実Dは既約被約quintic、S沿い位数は二とする。
最初のnormal quadratic G2∈H0(O_S(3))⊗Sym2 C2はgeneric squarefreeと仮定する。
generic double quadratic、m3、非plane等は別残余である。
本ノートに新しいJacobi/minor/GB/jet計算は無い。

## 1. discriminant-squareから全bihomogeneous factorizationへ

`PLANE_REDUCIBLE_CUBIC_CONNECTION_ENTRANCE.md` の全Hamiltonian parallel補題により、
degree二twisted normal representationのdiscriminantはinvariantであり、

\[
 \Delta(G_2)=c f^2,\qquad c\ne0.
\]

generic squarefreeの非零条件を明記する。Δ≡0枝をこれで捨てない。
K=C(S)上でbinary quadraticがsplitし、UFD C[x,y,z]のGauss lemmaにより

\[
 G_2=h\,\ell m,\qquad
 \ell=a u+bv,\quad m=c_1u+dv
\]

とfactor出来る。各row(a,b),(c1,d)はprimitive、coefficientsはそれぞれ同じbase degree d1,d2。
hはbase polynomial。bihomogeneous factorizationはbase/normalの二gradingを保つ。
rowの共通因子は全てhへ戻す。
二rowsはgenerically independentで、非零constantを吸収すると

\[
 h\det\begin{pmatrix}a&b\\c_1&d\end{pmatrix}=f. \tag{1}
\]

この式が全sourceであり、split rootをbasepointfree mapだと仮定しない。
rowの共通零点やR交点のcodimension二境界を保持する。
root交換でd1≤d2とすれば

\[
 (d_1,d_2)=(0,0),(0,1),(0,2),(0,3),(1,1),(1,2)
\]

の六degree flagsだけである。
f reducedの場合、hとdetは互いに素なR-component productsで、unused componentsを許す。
f nonreducedでも(1)とdegree flagsは有効だが、coprimalityだけは使わない。

## 2. triangular flagsとlinear/quadratic matrix flag

d1=0なら定数normal GL2でℓ=uと出来、

\[
 G_2=h\,u(A_{d_2}u+B_{d_2}v),\quad hB_{d_2}=f,
 \quad\gcd(A_{d_2},B_{d_2})=1.
\]

Aは全homogeneous polynomialを保持する。gcdの失敗はrow因子をhへ戻した別flagであり、
無断で補集合を消す条件ではない。
d2=0はconstant independent rootsを持つG2=κfuvへ戻る。
(1,1)はdetがf/hとなる全2×2 linear matrix sourceであり、
root rowsをmixするleft GL2は一般に積ℓmを変えるため、source parameterを保たずに除いてはいけない。

(1,2)ではhはconstant。primitive linear rowはindependent linear forms L1,L2であり、
point p=(L1=L2=0)は必ずR上にある。
逆にp∈Rならf∈(L1,L2)なのでL1D0−L2C0=fのquadratics C0,D0が存在し、全solutionsは

\[
 (C,D)=(C_0,D_0)+H(L_1,L_2),\qquad H\in H^0(O_S(1)).
\]

全H三係数とprimitive row条件の補集合を保持する。
Rがreduced reducible四型ならpのAut(R)orbitは有限である。
triangle/concurrentはsmooth-component点とintersection点の二型、
conic+line二型はconic上open、line上open、intersectionの三型である。
ただしこのfinite leading reductionからfull Poisson次元や帰属はまだ従わない。

## 3. triangle・constant-root flagは全高次sourceまでMPS

Rがnonconcurrent triangle f=xyz、degree flag(0,0)とする。
全normal/sourceを保つ定数GLでG2=κxyzuv、κ≠0。
既受理のglobal normal filtrationの全tensorをP=P0+P1+P2と書く。
P0, P1, P2はnormal weights0,1,2であり、

\[
 P_t=P_0+tP_1+t^2P_2
\]

はnormal座標の実GL5 rescalingから得る全Poisson曲線である。
t≠0では元のtensorとPGL同値、全coefficientsが多項式なのでt=0もPoisson。
weight零P0の実退化因子は正確にG2、従って

\[
 D(P_0)=\kappa xyzuv
\]

という五つのindependent coordinate hyperplanesのSNC divisorである。
κ≠0によりrank四log-symplectic SNC pointであることも確認出来る。

このirreducible curveの像全体を含む周囲Poisson成分Cを取るとCはSNC点P0を含む。
既受理のMPS Theorem1.4、すなわち **SNC log-symplectic locusと交わる全成分がMPS40**
という全component scopeにより、Cは既知MPS成分であり、元のP1点もその閉包へ入る。
適用scopeは `HYPERPLANE_QUARTIC_SEMISIMPLE_SUPPORT.md` §6、
`DOUBLE_QUADRIC_SUPPORT.md` と `KNOWN_COMPONENT_AUDIT.md` の一次文献監査と同じである。
SNC点が既知閉和に属するだけという誤った推論も、T01の一意性も用いない。

よって **triangle、m2、constant distinct normal rootsの全sourceはMPSに被覆**される。
normal connectionのresidueや全P1/P2を固定したsliceではなく、全sourceへのactual curveである。

## 4. 未完了の補集合

残る五degree flags、conic+line/concurrentのconstant-root endpoint、generic double normal quadratic、
m3、nonreducedRの残sourceは未完了。
nonreducedfやconcurrentfのendpointをSNCと誤認しない。
factorization (1)やR上pointの有限化だけを全source dimension/E帰属と呼ばない。
新しいfull-source検査は、必要最弱のwholeflag義務を選び別契約をgateしてから行う。
本ノートの新幾何と指定triangle全曲線は別作者独立読解を要する。
