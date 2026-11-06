# rank三の二次曲面：A=0の七cofactor型を六つの全磁気源へ還元

2026-11-07。課題 `E_reduced_nonnormal`。
対象はqがrank三、三次因子Cが既約であるD=qCの全支持のうち、
`DEGENERATE_QUADRIC_CUBIC_COFACTOR_REDUCTION.md` のA=0の七純冪零型。
全150線形源は共通の `CLASSIFICATION_DEGENERATE_QUADRIC_SOURCES_REVIEW_20261107.json` を入力とする。
本稿は七cofactor型の全支持を、最後の六つの全R6源へ帰着する。
§7–8の全周囲二次・真正三次障害により、全六源が相異なる23次元成分となる。
二次・三次障害を使う全六源の独立算術監査も完了した。
成分台帳はrootが別途更新する。

## 1. Cnil2の三型は全て空

### radicalへ非零の追加がある二型

全線形源のEuler収縮Pfaffianは、isotropic型ではx0q、nonisotropic型ではx1qで既に割れる。
従ってJacobiを課した実退化式の三次因子も可約であり、対象支持は空。

### 追加が零の型

U=(a,b,c)=(x0,x1,x2)、r=x3、s=x4、q=2ac+b²と置く。
三次元の線形bivector JをJab=a、Jac=−b、Jbc=cで定める。
全二十一源のU,r上のblockは

    R=(A+Br)J+rW∧∂r,
    A=−2n6a+2n9b−2n1c, B=n19, W=−J#dA/2.            (1)

実D=qCはunimodular liftのCasimirなのでCはRのCasimirである。
A≠0ならr成分からW(C)=0。
qのso3 Hamiltonian場Wの多項式不変環は、Aが非等方なら半単純場のweight計算、
等方ならJordan3の局所冪零場の不変環計算により、いずれもC[q,A]である。
従ってC∈C[q,A,r]。残るCasimir条件は

    [−2(A+Br)∂A+r∂r]C=0.

A'=A+2Br/3とするとA',rのweightは−2,1。
通常の同次次数が三でweight零の単項式はA'r²だけなので、Cは可約。
A=0,B≠0ならC∈C[q,r]の次数三部分はr³とqrで尽くされ、やはり可約。
A=B=0ならR=0、全五次元bivectorは一方向∂sとの外積でrank≤2。
この全境界を含む不変環論法は別担当も独立に受理した。

## 2. Czero、rank一非等方の型も空

同じqでcofactorはb∂rとなる。X不変性から残る四変数U,sのblock Rはr非依存で、bはR-Poisson-normal。
H:b=0上のq_H=2acは非平方である。
Cは既約三次なのでC_H≠0であり、q_HとC_Hは代数的に独立。
実際、同次二次・三次の代数従属はC_H²=κq_H³を意味するが、非平方q_Hに矛盾する。
H上の三次元skewテンソルが二つの独立Casimirを持つためR|H=0。
全源のliteral係数からn11=n17=n7=0が必要で、残るblockは

    R=n6 b(−a∂a+c∂c)∧∂s.

n6≠0ならC∈C[ac,b]なので同次三次式はbで割れる。
n6=0ならR=0となり全rank≤2。従って既約C支持は空。

## 3. Czero、rank一等方の型：fiber依存枝は空

cofactorはa∂r。全二十二源のU,s blockは、α=−n7/2と置いて

    R_UU=αaJ,
    a R_i,s=−Σ_j J_ij ∂j F,
    F=αa²s−n11a²b−n6a³/3+n15a²c−n18(abc+b³/3).       (2)

この全係数恒等式はJacobiを課す前から成り立つ。

α≠0の場合、Rの有理Casimir体はC(q,F)である。
aをinvertするとsとcをq,F,a,bから有理に復元でき、一般fiberは幾何学的に既約な有理曲面となるからである。
さらに写像(q,F)はa=1,b=0と置いてc=q/2、sをFから解く多項式sectionを持つ。
従って多項式Casimir環はちょうどC[q,F]。
実三次CはFの非零定数倍でなければならない。
n18=0ならFはa²で割れるので、既約C支持ではn18≠0。

sへのUの線形shearでF中のa²倍の三項を消し、
(a,b,c,r,s)→(λa,b,λ^(−1)c,νr,μs)と全bivectorの非零scalar倍を使えば、
q-cofactorをa∂rに戻したままα=1、F=a²s+b³+3abcへ規格化できる。
具体的にはshear後F/α=a²s+β(b³+3abc)、β=−n18/(3α)≠0であり、
ν=α、μ=(βλ²)^(−1)、scalar λ/νとすればよい。
従って全源で

    n7=−2, n18=−3, n11=n6=n15=0

を固定し、残る十七係数を全て自由に保持してよい。
この全源に三次Casimir式

    d(a²s+b³+3abc)·Π=−(a²s+b³+3abc)(a∂r)

を課す21×17の線形系はrank13、拡大rank14で不整合である。
`CLASSIFICATION_DEGENERATE_QUADRIC_Q3_AZERO_FIBER_20261107.json` に全入口を記録した。
後記の独立Fraction証明書は全二十一行と、それらの線形結合が1となる疎な双対証拠を保存する。
従ってfiber依存枝全体が空である。

α=0ではR=Y(U)∧∂s。rank四ならY≠0なのでC_s=0、従ってC=C(U)。
残るmixed列はU_0(U)+sL(U)で、Lは線形場。
q-normality及びD-CasimirからLq=LC=0。
L≠0なら前節と同じso3不変環からC∈C[q,linear]で可約となるためL=0。
従ってこの補集合も、次節の全plane-cone源に含まれる。

## 4. 残る全源はfree plane conic+cubicである

rank二のBを持つ二型は、全線形源でP_UU=0、mixed二列U,VはU変数だけに依存し、D=qCもUだけに依存する。
前節のrank一等方の補集合も同じ形に入る。
Dがrank四の実退化式なので

    det(E_U,U,V)=κD, κ≠0,
    U(D)=V(D)=0.

E_U,U,VはDのlogarithmic derivationの自由基底となる。
これはCramerで直接示せる。任意のlogarithmic derivationと一列を交換した行列式は、
Dの各既約因子上の一般点で零になるためsquarefree Dで割れ、Cramer係数が多項式になる。
従ってDは自由平面五次曲線でexponentsは(1,2,2)。
Jacobian idealのHilbert numerator 1−3t⁴+2t⁶から全Tjurina数はτ=12。

二次gradient syzygy空間は厳密に二次元で、任意のmixed基底はGL2のfiber座標変更で交換できる。
[U,V]=(div V)U−(div U)Vであるから、全Poisson源は

    Π=U∧∂r+V∧∂s+
      [−(div V)r+(div U)s+R2(a,b,c)]∂r∧∂s,          (3)

R2が任意の二次式、すなわち全六自由係数を保持した形で尽くされる。
この式はunimodularでもあり、その実Pfaffianは非零定数倍のqC。

## 5. 平面配置は六つに尽くされる

qは滑らかなconic、Cは既約cubicである。
kをqとCの異なる交点数とするとBézoutとMilnor数の積公式から

    Σμ(qC)=μ(C)+12−k.

τ=12≤Σμより、smooth Cは不可能。
Cがnodeならk≤1、cuspならk≤2である。

### cusp

C=b²c−a³、正規化(a,b,c)=(r²,r³,1)とする。
conicのpullbackは次数六でr¹項を持たず、逆にこの条件を満たす各多項式からconicが一意に復元される。

cusp r=0を避ける場合、交点が∞と一つの有限点ならr¹係数の条件に反する。
唯一の交点が∞の場合はpullbackが定数でq=c²となり、rank一なので除外される。
二つの非零有限交点を1とρ、重複度mと6−mにすると
ρ=−(6−m)/m。交換を除いてm=1,2,3の三候補。
m=2のconicは(3a+b−4c)²でrank一なので除外される。
一つのsmooth交点だけに重複度六で交わる候補はr¹係数が零にならず存在しない。

cuspを通る場合、その交点重複度はrank三conicなら2又は3である。
重複度四以上ではconicのcを含む全係数が零になりrank≤2となる。
もう一つの交点が∞ならpullbackはr²又はr³、q=ac又はbcとなりrank≤2なので除外。
有限のもう一つの交点をGm作用で1へ移すとpullbackはr²(r−1)⁴又はr³(r−1)³。
従ってcusp側に四つのrank三候補が残る。

### node

C=b²c−a²(a+c)、正規化(a,b,c)=(r²−1,r(r²−1),1)とする。
conicのpullback fはf(1)=f(−1)を満たす次数六多項式であり、逆にこの条件で一意に復元される。
唯一の交点が∞ならq=c²なので除外する。
唯一の交点が有限smooth点r=tならf=(r−t)⁶であり、
ζ=(t−1)/(t+1)は六乗根となる。
ζ³=1の三点はflex接線の二乗を与えrank一なので除外される。
ζ³=−1の三点は、nodeのprojective安定群ζ→ωζ（ω³=1）で交換され、t=0一代表でよい。
この作用がprojective線形なのは正規化の三次section空間が
span(ζ,ζ²,1−ζ³)であることから分かる。

唯一の交点がnodeなら、滑らかなconicは二枝の片方にしか接せないので、
二枝上の重複度は(1,5)又は(5,1)。nodeの枝交換で一代表となり、
f=(r−1)(r+1)⁵を採る。
従ってnode側は二候補。

### 全六源の具体的代表

| 名前 | q | C |
|---|---|---|
| cusp outside 1 | 225a²+24ab+1875ac+b²+1000bc−3125c² | b²c−a³ |
| cusp outside 3 | −3a²+3ac+b²−c² | b²c−a³ |
| cusp at cusp 2 | 6a²−4ab+ac+b²−4bc | b²c−a³ |
| cusp at cusp 3 | 3a²−3ab+b²−bc | b²c−a³ |
| node smooth contact | 2a²+3ac+b²+c² | b²c−a²(a+c) |
| node at node | 7a²+4ab+8ac+b²+8bc | b²c−a²(a+c) |

各代表で全十八係数の二次syzygy系を解き、kernel次元が厳密に二、
Saito determinantがqCの非零定数倍であることを確認した。
従って必要条件だけで作った候補から偽の源を加えてはいない。
全六源のq-cofactor Bのrankは二なので、§3のrank一等方の補集合も結局空である。

## 6. 証拠と残る義務

`CLASSIFICATION_DEGENERATE_QUADRIC_Q3_AZERO_REDUCTION_REVIEW_20261107.json`
はproducer/SymPyをimportしないGaussian Fraction実装で、全七源のPfaffian、
§1–3の全block恒等式、及びfiber-cubicの疎な不整合証拠を検査した（0.233秒）。

`CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_SOURCES_20261107.json`
は全七候補のsyzygy系、除外するrank一conic、残六候補の完全なanchor基底、
全六磁気係数でのJacobi・divergence・実退化式を保存する（1.60秒）。
全六配置の幾何尽くしは、∞に関する全境界も明示した後、rootが独立読解して受理した。

## 7. 全六磁気源のgeneric軌道と真の横断支持

(3)に対するfiber shear r'=r+L(U), s'=s+M(U)は、磁気項を

    R2 → R2+(V+div V)L−(U+div U)M

と変える。
全六源でこの6×6写像の階数は五であり、一つの商座標μが残る。
`CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_CANONICAL_20261107.json`
は全六磁気係数のまま、μK(U)への普遍shear行列を保存・検証する。
fiberをλ倍したpushforwardに全bivectorのλ^(−1)倍を合わせるとμはλμに変わる。
従ってμ≠0は各型で一つのPGL軌道、μ=0はその実族の極限である。
全R6源の飽和閉包は、このgeneric軌道閉包と一致する。

`CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_COMPONENTS_20261107.json`
は全150係数で、全六源のgroup・vertical・scalarの階数が48であり、
六磁気方向を全て加えても48であることを確認する。
従ってgeneric射影軌道次元は23である。
全Jacobi微分Dの階数と、Dの全kernelから真にgroup接を除いたnormalの二次支持は次の通り。

| 源 | rank D | 真normal数 | 二次余核の行空間 |
|---|---:|---:|---|
| cusp outside 1 | 101 | 1 | n0² |
| cusp outside 3 | 100 | 2 | n0², n0n1 |
| cusp at cusp 2 | 101 | 1 | 0 |
| cusp at cusp 3 | 100 | 2 | n0n1, n1² |
| node smooth contact | 101 | 1 | n0² |
| node at node | 101 | 1 | n0² |

二次の記録は `CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_NORMALS_20261107.json`。
各literal normal、全Dとの積零、Dの全150列を消す疎functional、
二次Jacobiとの非零pairingを保持している。
表の正次元支持を実変形と同一視しない。

一法線の非零二次障害を持つ三源では、実48方向に横断する48本の真の線形cutを取る。
残る102次元sliceの101個の一次式を陰関数で消すと、一変数完備環となる。
λD=0かつλJac(R)≠0は消去後も非零二次初項を与えるので、被約sliceは原点のみ。
48cutのKrull不等式から元の全150係数上の射影Jacobi局所次元は48以下であり、
実軌道48方向の下界と一致する。
従って各既約軌道閉包は全空間の23次元成分である。
この形式支持論理はgeneric_structureが独立に受理した。
局所唯一性、被約性、正確なnilpotence指数は必要なく主張もしない。

`CLASSIFICATION_Q3_FREE_CONIC_CUBIC_COMPONENT_REVIEW_20261107.json` は独立Fraction実装で
三型の全D101、G48、kernel49、101minor、全150列へのλD=0、非零二次pairingを再構成した。
さらに全六磁気係数のJacobi・divergence、全源のPGL shear、rank五の磁気商、
canonical尺度、μ=0の実極限も確認した（8.10秒、PASS）。

二次でfree軸が残る三型も、最初に可能な真正三次障害によって閉じる。
free軸のnormalをRとし、全Dの像でDT=−Jac(R)を解く。
三次項B(R,T)は単独で使わず、全normalの全二次pairing Qijを加えた空間

    im D + span(Qij)

で割る。
`CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_THIRD_20261107.json`
は三型全てでλD=0、λQij=0、λB(R,T)≠0を示す。
さらに全51個のactual group・Euler-vertical・scalar列G_jについて
B(R,G_j)が同じ空間に入ることも厳密に検査した。
従ってTを任意のfull kernel方向で変えても、この余核の非零値は変わらない。
特定の二次liftの残余だけを障害と呼んでいるのではない。

全normal二次を消すλに対応するIFT後の残式は次数二以下が零であり、
free軸上の三次係数が非零となる。
outside3、atcusp3では初期二次idealの零点が一軸であり、
その軸上非零の三次式を加えると接錐支持は原点だけになる。
atcusp2は一変数で最初の三次式が非零なので同じ結論となる。
前記と同じ真48cut/Krullの議論から、残三つの軌道閉包も23次元成分である。
この真正三次の形式論理はrootとgeneric_structureが独立に受理した。
`CLASSIFICATION_Q3_FREE_CONIC_CUBIC_THIRD_REVIEW_20261107.json` は独立Fraction実装で
残三型の全源、全D/G、真正二次行空間、DT+Jac(R)=0、λD=λQij=0、
全51個のB(R,G_j)の消滅を再構成した（8.32秒、PASS）。
非零三次pairingはoutside3で−4480/9、atcusp2で547805211540848640000000/211、
atcusp3で−241920/11となる。

途中で指定したP+tR+t²TがそれだけではPoissonでなかった診断
`CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_CURVES_20261107.json`
は、Tを保存する入口としてのみ用いる。
高次jetからの外挿、実変形の不存在、正確なnilpotence指数は主張しない。

## 8. 六つの23次元成分は互いに、かつ既知74成分全てと異なる

六つのgeneric軌道は以下の幾何型を持つ。

* cusp cubic、conicはcuspを避け、二つのsmooth交点に重複度1と5で交わる。
* cusp cubic、conicはcuspを避け、二つのsmooth交点に重複度3と3で交わる。
* cusp cubic、conicはcuspで重複度2、もう一つのsmooth交点で重複度4。
* cusp cubic、conicはcuspで重複度3、もう一つのsmooth交点で重複度3。
* node cubic、conicはnodeを避け、唯一のsmooth交点に重複度6で交わる。
* node cubic、conicはnodeだけで交わり、その二枝の重複度は1と5。

実五次式の本質三変数空間はintrinsicに復元できる。
従ってPGL5による同値は平面上のqとCの配置も同値にし、上の六配置は互いに同値でない。
各成分は23次元の稠密軌道を持つので、二つの成分が一致すれば二つの稠密軌道が開集合で交わる。
それは不可能であり、六成分は相異なる。

全六閉包のrank四部分では、退化五次式の本質変数数が三以下という閉条件が成り立つ。
この条件は五つの一次偏微分の係数行列のrank≤3という行列式条件である。
全MPS40成分は五つの独立な超平面の積を退化式に持つSNC点を含み、その本質変数数は五。
従って全MPS40成分のいずれとも一致しない。
単にgeneric点が非SNCであることを区別の根拠にはしていない。

六rank二成分とは一般Poisson rank、他次元の既知rank四成分とは次元で異なる。
残る既知23次元非MPS成分の一般因子型は
`KNOWN_COMPONENT_AUDIT.md` と `ADDITIONAL_RANK4_LATE_COMPONENT_AUDIT.md` の表で尽くされる。
H_QUARTIC_ESSENTIAL5、MINIMAL_M、A0、DPMはH＋既約四次、
ZERO_LOG_5Hは5H、CMINUS_MAG及び独立三H・三接線二次は三一次因子＋既約二次である。
いずれも本稿のgeneric既約二次＋既約三次と異なる。
同じ成分なら両generic constructible像の非空開が交わり因子型が一致するはずなので矛盾する。
並行研究の二一次因子＋既約三次型の新23成分候補とも、同じgeneric因子型の議論で異なる。

以上は六成分を成分台帳へ接続するための相違証明である。
全六成分の独立算術監査は完了しており、A=0の七型の全支持はこの六成分に被覆される。
全分類上界は主張しない。
