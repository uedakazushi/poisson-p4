# 非退化 member を持つ integral CI22：独立幾何監査

2026-11-07、主義務 G。**PASS**。
`REGULAR_PENCIL_CI22_ABELIAN_CONORMAL_NONDOMINANCE.md` の強い最終版を確認した。
対象は integral・非退化な S=V(q1,q2)⊂P4、pencil に非退化 member があり、
π|S=0、実退化因子 D が既約被約 quintic である全 rank4 Poisson tensor。
全 PGL 飽和像の次元は≤24。S の smoothness、normality、node、clean 条件を使わない。
全 member が退化する pencil は残す。全 E 帰属又は rank4 不存在は結論しない。

## 1. 二十源と abelian conormal 条件

q1,q2 が regular sequence なので、E=Λ²TP4 について
`0→E(-4)→E(-2)²→E I_S→0` は特異 S でも exact。
外積 Euler と line-bundle cohomology から H0(E(-4))=H1(E(-4))=0、
H0(E(-2))=Λ²C5 であり、全源は一意に Π=q1A+q2B、A,B constant skew matrix。

π|S=0 なら D は S で少なくとも二重。三重以上なら一般二点を結ぶ secant line は
次数五の D に全て含まれる。非退化非plane 曲面の secant variety は P4 全体なので矛盾。
この secant 事実は、dimSec≤3 とすると一般 p の三次元 join(p,S) が Sec と等しくなり、
一般 p が全て Sec の vertex になることからも従う。vertex locus は線形で、S が P4 を張るので矛盾する。
従って generic normal multiplicity は正確に二。

一般滑らかな点で π の linearization は、二次元 conormal ideal N* と abelian 商 C² を持つ
四次元 Lie algebra。非零 normal quadratic Pfaffian はこれを Frobenius にする。
N* が aff1 なら Der(aff1)=ad(aff1)、center=0 より商の二 lift を調整して
aff1 と可換に出来、両 lift も可換になる。従って aff1⊕C²、generic rank≤2 となり矛盾。
よって N* は abelian で、generic に {I,I}⊂I²。
この議論は `RANK_ZERO_SURFACE_GENERIC_FIRST_JET.md` §§1–3 とも独立に照合した。

prime complete-intersection I=(q1,q2) では gr_I R=(R/I)[t1,t2] は domain。
a∉I と b∉I² の最初の I-adic symbol の積は零にならないため I² は I-primary。
従って generic symbolic-square 条件は全 homogeneous ordinary square 条件になる。
projective 局所式 qi/ell² への補正も、Π の各係数が I にあるため I² に入る。
{q1,q2}=q1hA+q2hB∈I² と I/I² の自由性から hA,hB∈I2=span(q1,q2)。
ここは必要線形条件であり、全 Jacobi と同値とは主張しない。

## 2. normalizer と全 Jordan 型

q1=xᵀQx、q2=xᵀQMx、Q 非退化 symmetric、M self-adjoint とする。
T=AQ は Q-skew で
`dq1 A dq2 = 2xᵀQ[T,M]x`。従って A,B はともに
`N(M)={T∈so(Q):[T,M]∈span(I,M)}` に対応する。

異なる固有値が二つ以上なら、各 generalized eigenspace への trace は
α+βλ=0 を与え、[T,M]=αI+βM は零となる。
一固有値なら scalar を引いて N nilpotent とし、trace から α=0。
中心化した Jordan grading S が [S,N]=N を満たすので、normalizer は skew centralizer⊕CS。

どの λ でも q2−λq1 は rank≥3。rank≤2 なら複素数体上で一次式の積になり、
I の素性から S が一超平面に含まれてしまう。
従って各 λ の geometric multiplicity は≤2。

self-adjoint Jordan 型は Q-orthogonal cyclic block に分かれる。
これは最大 nilpotent 長 m に対し Q(v,N^(m−1)v)≠0 となる v を選び、
その非退化 cyclic subspace と直交補空間へ分ける帰納法で得られる。
各 block の Q は unit polynomial の平方根による generator 調整で anti-diagonal に出来る。
複素数体上で余分な符号不変量はない。
一 block の commutant は C[N] で self-adjoint なので skew 部分は零。
二 block 長 s,t の間は min(s,t) 次元の intertwiner とその skew transpose のみ。
ゆえに skew centralizer の次元 c=Σmin(s,t)≤2。

## 3. 二固有値以上を全て排除する common image

同一 λ の二 block 長 s≥t の間の全 intertwiner は
Wλ=ker(Nbig^t)⊕small block に像を持つ。dimWλ=2t。
他の固有値間の項と一 block 内の skew 項は零なので、全 centralizer の像は
固定 W=⊕Wλ、dimW=2c≤4 に入る。A=TQ⁻¹ の像も同じ W。
従って全 Π の像は W、すなわち Π∈Λ²W。
dimW≤3 なら Π∧Π=0。dimW=4 なら Λ4W は固定一次元であり、
実 D=Euler∧Π∧Π は W を消す固定一次式 ell(x) で割れる。
これは全 eigenvalue/crossratio/Jordan degeneracy に対する構造証明であり、標本計算ではない。

## 4. 一固有値の三型と残る六源

geometric multiplicity≤2 より partition は 5、4+1、3+2 の三つ。
J5 は normalizer 一次元なので Π=α(q1,q2)C、D は α² と一次式の積又は零。
J4+1 は normalizer が grading S と centralizer generator R の二次元。
rankR=2 より bivector C=RQ⁻¹ は C∧C=0。
Π=α SQ⁻¹+β C の外積平方は α で割れ、D も非零 quadric α を因子に持つ。
**matrix R²=0 は使わない。実際に後述の canonical matrix は R²≠0 である。**

唯一残る J3+2 は normalizer 三次元、全必要線形源は六次元。
R を5×5 superdiagonal shift、Q を anti-diagonal (+1,−1,+1,−1,+1)、
N=R²、S=diag(1,1/2,0,−1/2,−1) とする。
R,R³ は skew centralizer の basis、S は skew grading。
K=(NS+SN)/2、L=N/2+K に対し

`δ_R(q1,q2)=δ_R³(q1,q2)=0`,
`δ_S(q1,q2)=(0,−q2)`, `δ_L(q1,q2)=(q2,0)`。

従って四つの線形 field は actual pencil-preserving algebraic group の Lie algebra に入る。
char0 の閉 algebraic group は smooth なのでその次元は少なくとも四。
各 field の flow としても、R,R³,L は nilpotent exponential、2S は整数重み C* に積分する。
これらの線形結合が scalar I なら Q-self-adjoint part から L の係数と scalar が零、
N との commutator から S の係数が零、R,R³ の独立性から残りも零。
未証明の全 Möbius 群の lift は不要である。

六源上の Poisson 条件は proper 閉条件である。
実際 A=SQ⁻¹ は rank4 で Π=q1A は源に入るが projective Poisson でない。
一般に rank4 constant A の kernel covector を x4 とし、残る四次元を symplectic にする。
E∧[qA,qA]=0 なら x4 を含む成分から q A♯dq∧A=0。
四次元では v↦v∧A が同型なので A♯dq=0、q は kernel coordinate の平方でなければならない。
非退化 q1 に矛盾する。よって全 Poisson 源の各成分は次元≤5。

J3+2 pencil は上の self-adjoint 標準形により一つの PGL orbit。
必要条件と Poisson 条件は intrinsic なので実 preserving group H（dim≥4）が全源 X を保つ。
PGL5×(X\{0})→P H0(E) の fiber には
`(h,λ)↦(g h⁻¹, λ h·Π)` という五次元部分が入る。
h は第一成分から決まり、Π≠0 より λ も決まるため、個別 Π の stabilizer を仮定しない。
従って image 次元≤24+5−4−1=24。
これで残余型全体の G 非優越が従う。

## 5. 有限 QQ 補助証明書

`code/classification_regular_pencil_ci22_review.py` は producer/CAS import なしの独立 Fraction 実装。
三つの全単固有値 normal form に対し、全10 skew A 係数と α,β の25×12行列
`[AQ,N]−αI−βN=0` を構成し、全 rank と非零 minor、完全 basis を検査した。
kernel 次元は順に1,2,3。J4+1の rankR=2、C∧C=0、matrixR²≠0 も検査した。
J3+2の四 field の全 pencil 微分、scalar を含む五行列の rank5、全六源の rank6、
q1SQ⁻¹ の全350 projective Jacobi と非零 witness を保存した。
これらは全 pencil の代わりに三点を調べた推論ではなく、上の全標準形証明を補う有限算術である。

現行証明書は `reports/CLASSIFICATION_REGULAR_PENCIL_CI22_REVIEW_V2_20261107.json`、
契約は `reports/CLASSIFICATION_REGULAR_PENCIL_CI22_REVIEW_JOB_20261107.json`。
最初の task ID 誤記で gate が拒否した後に orchestration の誤りで短い probe が走ったため、
その初版 report は未受理の探索記録として保存した。
正しい active G 契約の gate PASS **後**に全計算を再実行した V2 だけを今回の有限証拠とする。
0.0298秒、最大 RSS 12,672KiB、全項目 PASS。
入力ノートの確定 snapshot、SHA256、実装全文、全行列と minor を結果に保持する。

全 member が退化する Kronecker pencil、common radical を持つ cone は未処理補集合。
本結論はその補集合を閉じず、他の曲面又は分類全体の完成を主張しない。
