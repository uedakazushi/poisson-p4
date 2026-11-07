# m4 generic31：全実quinticを二つのessential-four正常形へ

2026-11-07。主義務 G。rank二planeに沿って実既約被約quintic D の位数が四、
leading binary quarticのgeneric partitionが3+1、coefficient image rankが二とする。
ここでは zero cubic R の被約性や接続のlog性を仮定しない。
全支持を二つの actual-D normal formsへ還元し、必要最弱の全linear Casimir sourceを指定する。
新しい source dimension/rank 計算はまだ行っていない。

## 1. 全 actual D の二正常形

`PLANE_M4_REPEATED_LEADING_GEOMETRY.md` §2 の P1 Bertini により、
generic triple rootは全linear quartic spaceの固定triple rootである。
定数normal GL2とindependent base linear forms A,B を選ぶと

\[
 D=u^3(Au+Bv)+g_5(u,v).
\]

実Dの既約性からg5(0,v)=κv5、κ≠0である。
κ=0ならg5もuで割れ、D自体がuを因子に持つためである。
A,B へ normal linear forms を加える可逆linear shearで、g5のu5,u4v,u3v²を除く。
残る項をαu²v³+βuv⁴+κv5とする。
v↦v−βu/(5κ)でuv4を消し、再びu-degree≥3の項を同じlinear shearsへ吸収出来る。
u²v³の新係数はα−2β²/(5κ)。
全体尺度とA,Bの共通尺度でκを一にし、uの非零尺度で残る係数を零又は一へ戻す。
必要な平方根はC上の有限選択であり、全sourceを減らさない。

従って全actualDは、GL5変換と非零全体尺度を除いて

\[
 F_\epsilon=u^4A+u^3vB+\epsilon u^2v^3+v^5,
 \qquad \epsilon\in\{0,1\}
\]

の二つに尽くされる。残る第五座標をwとする。
両FはAに関してprimitive linearなのでC上絶対既約。
四つのpartial derivativesはC上独立で、essential-variable rankは四である。
ここでは二正常形の相違、存在するPoisson構造、成分性又は既知親族帰属を主張しない。

## 2. F0は全rank四sourceが空

base variables A,B,u,v上のhomogeneous quadratic gradient syzygyを
δ=P∂A+Q∂B+R∂u+S∂v とする。
F0のgradient式は

\[
 u^4P+u^3vQ+(4Au^3+3Bu^2v)R+(Bu^3+5v^4)S=0.
\]

u=0からS=uL、Lはlinear。式をuで割り再びu=0とするとL=cu、cはconstant。
もう一回uで割ってu=0とすると3B R|u0+5cv³=0。
Bとvが互いに素なのでc=0、R=uM、Mはlinear。
残る式はu(P+4AM)+v(Q+3BM)=0であり、
互いに素なu,vからP+4AM=vN、Q+3BM=−uN、Nはlinear。
したがって全quadratic gradient moduleは厳密に

\[
 \delta=N Z+M Y,\quad
 Z=v\partial_A-u\partial_B,\quad
 Y=-4A\partial_A-3B\partial_B+u\partial_u
\]

で尽くされる。全linear gradient syzygiesはZ,Yのconstant span、degree zeroは零である。

canonical divergencezero strict quadratic lift Π のPfaffian-Casimir identityによりΠ#dF0=0。
各rowをwについて展開すると、base成分はdegree0/1/2 gradient syzygiesである。
ゆえにbaseへ投影したΠの全imageは、pointwiseにspan(Z,Y)内にある。
vertex方向∂wを戻してもimage次元は高々三、skew tensorのrankは偶数なので高々二。
従って **F0をactual nonzero quinticに持つrank四Poisson tensorは存在しない。**
Rやlog性、normal connectionの一意性、full nonlinear Jacobiは不要である。

## 3. F1に対する必要最弱の全linear source

F1に対し五変数の全150quadratic bivector coefficientsを保持し、

\[
 L_1=\{\Pi:\operatorname{div}\Pi=0,\ \Pi^\sharp dF_1=0\}
\]

というC-linear subspaceを取る。
ここでdiv0はcanonical strict liftの正当な条件であり、任意Euler代表をstabilizerと置いていない。
Casimir identityの既受理の証明は `ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md` §1。
そのidentity自体はFのessential-variable数に依存しない。
固定actualD class [F1]を持つ全projective Poisson sourceはP(L1)内に入る。
JacobiとactualPf条件を課すことはこの必要源を狭めるだけである。

Z=v∂A−u∂BはF1も保存する非自明unipotent PGL4作用を与える。
したがって[F1]のPGL4 orbitは次元高々14。
essential-four base選択Gr(4,5)の次元四を戻し、全actualDのfamilyは次元高々18。
canonical liftのGL5 equivarianceによりdimL1はこの全orbit上一定なので、
全projective Poisson像は

\[
 \dim\le18+(\dim L_1-1).
\]

特に **dimL1≤7なら全m4 generic31 image-rank二sourceのG像≤24**が従う。
この場合、全Jacobi/実Pfaffian方程式のnonlinear eliminationや八接続profileのrank計算は不要。
dimL1が大きければ上界だけを保存し、根拠なくG完了や新成分とは呼ばない。
whole-sourceの実normal family又は他の幾何を設計する前に相談する。

## 4. exact finite jobの入口と補集合

必要な一回のjobはε0/1の全gradient-linear sourceを全150coefficientsで作るだけでよい。
五変数degree6の全係数によるΠ#dF=0とdivの全linear係数を保持し、
rank minor、完全kernel、全row代入を保存する。
F0については上のAR2 basisとの全source包含を別に確認出来る。
F1についてdimL1≤7なら§3の全source boundへ直結する。
単一点のambient tangent rankやjetをこの代用にしない。

image-rank一、generic22/4は既essential-three全八成分へ還元済み。
generic211は別六型packet。m2/m3、zeroR、他supportを勝手に消さない。
この二normal form還元はgeneric31 rank二に限るが、**Rのnonreduced枝も含む**。
全必要sourceの幾何と有限rank結果は独立reviewを要する。
