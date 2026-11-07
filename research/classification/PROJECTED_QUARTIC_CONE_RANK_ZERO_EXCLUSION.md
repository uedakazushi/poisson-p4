# projected rational quartic curveのcone：全centre parameterでrank零・既約quinticを排除

2026-11-07、主義務G。
minimal-degree surface S(0,4)⊂P5はrational normal quartic C4⊂P4のconeである。
centreをvertex shiftでbase P4へ移した後、centre pがclosed Sec(C4)の外にある場合を扱う。
image S⊂P4は、pから同型射影したsmooth rational quartic C⊂P3のconeである。

**命題。** この全一parameter family及び全collision型を含め、S上で零となるrank四Poisson tensorの
実Dは既約被約quinticにならない。
generic centreを一つ選ぶ計算は用いず、nonCMなextra sectionも全て保持する。

独立reviewでの最弱の証明は、§2の全source表示から直ちに得られる。
すなわち `Π=Q*A+∂c∧V` なら

`Π∧Π=Q*(Q*A∧A+2*A∧∂c∧V)`。

従って全raw four-vector係数がQで割れ、Eulerとのwedgeから定めるactualDにもQ因子がある。
既約quintic排除にはJacobi、conormal条件、centre stabilizerの計算を必要としない。
§3–5は追加の強い帰結のproducer導出として残すが、G排除の証明・reviewは§1–2とこの因子恒等式だけで完結する。

## 1. base curveのunique smooth quadric

C≅P1、O_C(1)=O_P1(4)なので、h⁰(O_C(2))=9、ambient quadricは十次元。
従って少なくとも一つquadric QがCを含む。
二つ独立ならCはdegree四CI(2,2)でarithmetic genus一となり、smooth rational Cに反する。
よってQは一意である。
rank二以下ならCのintegralityと非退化性に反する。
rank三quadric coneのclass groupはruling generatorRについてH=2Rであり、degree四Cのclassは4R=2H。
従ってこの場合もCはquadric sectionでCI(2,2)となるので不可能。
よってQはrank四のsmooth quadric surfaceである。

## 2. vertex weightで全projective sectionを尽くす

base座標をx=(x0,x1,x2,x3)、vertex座標をcとする。
全global bivectorのraw homogeneous二次liftについて、cのweightを一、baseを零とする。
一般形はbase-base部分 `P2(x)+c L1(x)+c²A0` と、
base-vertex部分 `∂c∧(V2(x)+c W1(x)+c²v0)` である。
全cone上でprojective sectionが零という条件は
`(Euler_base+c∂c)∧Π(x,c)=0` を全Cと全cに課すことである。

weight二のc³係数からA0=0。
weight一のc²係数から、L1はC上でEuler_base∧constant vectorと等しくなる。
CがP3をspanするためこのlinear-coefficient equalityはambient全体に成立し、
weight一部分はEuler_full∧linear fieldとして全てgaugeで除ける。

weight零では、P2からEuler_base∧W1を引いた全quadratic係数がC上で零となる。
I_C(2)=C*Qなので、この部分はQ*A、A任意constant skew4×4 matrixである。
対応するmixed c項も同じEuler_full gaugeで同時に除ける。

従って**全projective rank零section**は、適切なEuler liftで

`Π=Q(x)*A+∂c∧V(x)`                                     (1)

となる。Vは任意のquadratic vector fieldで、C上ではEuler_baseとproportionalとなる条件だけを持つ。
VをQ*constant vectorに制限していない。
特にCがP3でlinearly nonnormalであることから生じるextra sectionも、この全radial-on-C Vに含まれる。
vertex nonnormality又はh¹(O_S)を無視した十係数ansatzではない。

## 3. full projective Jacobiはstrict base Jacobiを強制する

(1)はcに依存しない。strict bracketを取ると

`[Π,Π]=[QA,QA]+2∂c∧[V,QA]`。

Euler_fullとのwedgeのc係数は `c∂c∧[QA,QA]` だけなので、
全projective Jacobiは `[QA,QA]=0` を強制する。
もしrankA=4なら、四次元symplectic線形代数で

`[QA,QA]` は非零定数を除いて `Q*A^sharp(dQ)∧A`。

v→v∧Aは同型なのでA^sharp(dQ)=0となるが、Qは非退化であり不可能。
従ってPoisson点ではrankA≤2である。

## 4. Frobenius conormal条件はAのlinear fieldをCのautomorphismへ制限する

Sはnonplaneなので、既約quinticを仮定すればsecant/Frobenius補題からgeneric conormal idealはabelian。
Cのgeneric pointで、Qともう一つのhomogeneous polynomial Fがnormal basisを作るようにFを選ぶ。
global generator degreeを三に限定する必要は無い。degreeF=eとする。
vertex chart c≠0でQ/c²,F/c^eを使うと、projective bracket numeratorは

`c{Q,F}−eF{Q,c}−2Q{c,F}`。

(1)のbase部分からの寄与は `c*Q*(dQ A dF)`、vertex部分の寄与はcに依存しない。
conormal basis Q,Fに関してc係数を比較すると

`(dQ A dF)|C=0`。

従ってlinear field `X=A^sharp(dQ)` はQとFの双方にtangentであり、Cのgeneric pointにtangentとなる。
polynomial idealのintegralityにより全Cへtangencyが延び、Xはembedded Cのinfinitesimal projective automorphismである。
matrixとしてXはQ-skewであり、scalar Euler成分は持たない。

Vの寄与を捨てたのではない。
VはC上でradialなので、V(Q),V(F)はともにI_Cへ入り、vertex側の二項は既にI_C²へ入る。
このこともconormal比較と整合する。

## 5. 全binary quartic centreのstabilizerを処理する

embedded CのautomorphismはnormalizationP1へliftし、
四次元linear subsystem `W⊂H⁰(P1,O(4))` を保存する。
これはprojection centreのdual line pを固定するPGL2と同じ群である。
Sym⁴表現のself-dualityにより、そのconnected stabilizerはpをbinary quarticとして根で分類出来る。

根のsupportが三点以上なら、connected PGL2 stabilizerは各点を固定するため零次元。
この全caseではX=0、従ってA=0である。
根が一点ならpはC4自身に入り、centreとして除外される。
根が二点ならmultiplicityは3+1又は2+2。
3+1型はC4のtangent lineに入り、closed Sec(C4)内なので今回の仮定から除外される。

唯一残る全caseは2+2型で、一つのPGL2 orbitである。
centre pをC4のmiddle coordinateへ取れば

`C=(s⁴,s³t,s*t³,t⁴)⊂P3`, `Q=x0x3−x1x2`。

connected stabilizerはdiagonal scalingの一次元で、Q-skew representativeは

`Z=diag(2,1,−1,−2)`。

そのrankは四である。
従って、非零なQ-skew infinitesimal automorphism XはZのscalar倍であり、
`rankA=rankX=4` となる。これは§3のrankA≤2に反する。
よってこの全special orbitでもA=0。

全centre parameterでA=0となり、(1)はdecomposableな∂c∧Vだけとなる。
global rankは二以下で、rank四・既約quinticの仮定に矛盾する。命題が従う。

## 6. secant内centreとdegree四曲面への接続

centreがSec(C4)内の場合、coneのsecant locusは一dimensionalである。
直接読解済みのBrodmann–Parkのminimal-degree projection criterion
（`PROJECTED_QUARTIC_SCROLL_ISOLATED_MODELS.md` の一次資料参照、§3及び§6）により
imageはACMなdegree四codimension二surface、従ってCI22である。
その全rank零G枝は別証明で処理済みである。

本noteとsmooth S13/S22四model、Veroneseの二centre型、全CI22を合わせれば、
almost-minimal-degree分類の入力の下でdegree四rank零曲面の残projection branchesを全て処理出来る。
その全degree四分類への最終外部入力の適用範囲は独立reviewで別に確認する。
rank二restriction、degree五・六の非normal曲面、及びplaneはこの命題の補集合である。
