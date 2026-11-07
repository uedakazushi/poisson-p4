# normal CI(2,2)上のrank二：全tangent sourceと三つの有限残型

2026-11-07、主義務G、producer plane_collision_review。
本稿はnormal integral nondegenerate S=CI(2,2)と、実Dが既約被約quinticの全rank二restrictionを扱う。
non-normal CI22は明示的な別補集合である。

## 1. 全二十五tangent source

Sはnormal Gorenstein、ωS=OS(−1)。Sreg上のrestriction bivectorはOS(1)のsectionであり、
normal性とcodim二補集合により全Sへ延びる。従って全H0(OS(1))五係数を持つ。

U∈Λ4C5をconstant four-vectorとし、JU=ι(dq2)ι(dq1)Uを取る。
q1,q2はstrict Casimirであり、全Uは四vectorとしてdecomposableなのでJUはstrict Poisson、JU²0。
U=ι(ell)volumeと表すと、E∧JUのq1=q2=0 restrictionは非零constant倍の
ell(E)ι(dq2)ι(dq1)volumeとなる。
従って全五linear restrictionがJUによりliftされる。

差はS上零である。Koszul列とE=Λ²TP4に対するE(−4)=Ω²(1)のacyclic性から、
全差はq1A+q2B、A,B constant skewの全二十係数である。
従って全tangent tensorは

`Π=JU+q1A+q2B`。

JUは両qをstrict Casimirとするので、conormal-abelian必要条件は既rankzero sourceと同じで、
dq1 A dq2,dq1 B dq2∈span(q1,q2)となる。
同時tangent Euler liftが既に取れているので、local projective bracketの補正もI²へ入り、
prime CIのI²-primary性で全homogeneous条件になる。

## 2. regular pencilの全型を三つに縮める

非退化member q1=xTQxを持つとき、既受理の全self-adjoint Jordan補題を使う。
許されるAの全空間はN(M)、二固有値以上ではskew centralizer、単一固有値ではcentralizer＋grading。
各固有値のgeometric multiplicityは二以下、skew centralizer dimension c≤2。

k=dimN≤1では、Π=JU+alpha(q1,q2)Cとなる。
JU²0ゆえactualDはquadratic alphaで割れ、又は零。従って全既約quintic枝は空。

normal性から、各λのP(ker(M−λ))とq1の交わりは有限集合でなければならない。
二つのJordan blockの長さが双方二以上ならq1はその二次元kernel上恒等零となり、
singular lineを生じるのでnormalではない。
全dimension五とk≥2を組み合わせると、normalな全残型は次の三つである。

| model | q1 | q2 |
| --- | --- | --- |
| sem221 | x0²+x1²+x2²+x3²+x4² | x2²+x3²+2x4² |
| mixed21_11 | 2x0x1+x2²+x3²+x4² | x1²+x3²+x4² |
| nil41 | 2x0x3+2x1x2+x4² | 2x1x3+x2² |

全三固有値又は二固有値の位置はpencil PGL2で固定出来、単一固有値は零へ移せるので、
各々一つのactual PGL5 surface orbitである。
残りc=2の(22)+(1)と単一(3,2)は上のsingular-lineでnon-normal補集合に入る。
各表のNは二次元なので、全必要sourceはJU五係数＋q1,q2に二generatorの全四係数、計九。
後四係数の二×二matrixがrank一以下なら再びquadratic factorとなる。
rank二openだけが未処理であり、normalizer/元Jacobiから最弱のG上界又はactualD因子を狙う。
この時点で全sourceのPoisson支持を解いたとは主張しない。

## 3. 全member退化のnormal補集合はrank二以下

既受理のall-singular-pencil canonical coverを用いる。
M1+regular2の唯一integral modelはsingular lineを持つためnormalでない。
normalとして残るのはcommon radical一、smooth elliptic quartic C⊂P3のconeだけ。
base Cがsingularならそのconeはsingular lineを持つのでnormalでない。
C smoothならq1=Σxi²、q2=Σλixi²（i0..3、λ全相異）へ直せる。

cone vertex coordinateをvとする。conormal条件から全A/Bのbase partは零、
全sourceはJU+∂v∧(q1V+q2W)、V,W constant base vectors。
Uのbase-four-vector係数をcとし、J=Jbase(q1,q2)と置く。
c0なら全Πは∂v∧Yでrank二以下。
c非零では全UのcompatibilityからJUを含む元Jacobiは正確に

`Ebase∧(q1 Lie_V J+q2 Lie_W J)=0`

を必要とする。missing dのcomponentでxi⁴の係数は非零共通定数を除き

`(λd−λi)*(Vd+λi Wd)`、i≠d。

相異なる二iを取るとVd=Wd=0、全dについてV=W0となる。
従って全Π=JUで再びrank二以下。
この証明は全elliptic modulusと全vertex sourceを保持する。

## 4. 一回の有限whole-source診断

三表について全raw25/Euler rank、全N二generator、全九source/all150を復元し、
全350original Jacobiと実126Pfを保存する。
quadratic coefficient rowspaceをexact Fractionで取り、全equivalent rowspace basisをfactor表示する。
これはnormalform点の標本ではなく、上記三surface orbit上の全parameter sourceである。
GB、より高いjet、少数のfree-coefficient点による全支持推論はしない。
factored全元式からwhole rank二openの因子又はdimensionを手で閉じられなければ、
その全sourceを保存し、別の幾何経路へ変える。
