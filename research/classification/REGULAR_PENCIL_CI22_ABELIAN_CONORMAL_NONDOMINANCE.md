# 非退化memberを持つ全integral CI22：abelian conormal normalizerによるG上界

2026-11-07、主義務G。S⊂P4をintegral・非退化なCI(2,2)とし、
defining pencilに非退化quadricがあると仮定する。
**命題。** S上で零となり、実退化因子が既約被約quinticである全rank四Poisson tensorの
projective PGL5像はdimension二十四以下である。
Sのsmoothness、normality、node、clean条件を仮定せず、全self-adjoint Jordan型を扱う。
全memberが退化するpencilは明示的な補集合である。

## 1. 全二十係数とconormal bracket

Koszul/Euler証明により全rank零sourceは一意に
`Π=q1 A+q2 B`、A,Bは任意のconstant skew5×5 matrixで尽くされる。
これは特異CIでもregular sequenceだけで成立する。
nonplane Sと既約quinticという仮定から、secant/Frobenius補題により
Sのgeneric conormal idealはabelianである。

prime complete-intersection idealI=(q1,q2)のsquareはI-primaryである。
例えばgr_I R=(R/I)[t1,t2]がdomainなので、a∉Iなら
ab∈I²⇒b∈I²を得る。従ってsymbolic squareはordinary squareと一致する。
projective bracketのEuler補正はI²へ入るため

`{q1,q2}=q1 h_A+q2 h_B ∈I²`,

`h_A=dq1 A dq2, h_B=dq1 B dq2`

となる。I/I²がR/I上でq1,q2をbasisとするfree moduleなので、
h_A,h_BはそれぞれIに入る。
両者はdegree二だから、どちらもspan(q1,q2)に入る。
これは全sourceに対する**線形**必要条件であり、Poissonの全Jacobiをまだ減らしていない。

## 2. 条件はskew normalizerそのものとなる

`q1=x^T Q x`、`q2=x^T Q M x` と書く。Qは非退化symmetric、MはQ-self-adjoint。
Aに対してT=AQと置くとTはQ-skewであり、直接のmatrix計算から

`h_A=2*x^T Q[T,M]x`。

よって許される全Aは

`N(M)={T∈so(Q): [T,M]∈span(identity,M)}`                  (1)

と一対一に対応する。Bも独立に同じ空間を取る。

rをMの相異なる固有値の数とする。
r≥2なら、commutatorがαI+βMである場合、trace polynomial identitiesにより
各固有値λについてα+βλ=0となる。
従ってα=β=0で、N(M)はskew centralizerだけである。
r=1ではMからscalarを引いてnilpotentと出来、traceからα=0。
各Jordan blockの中心化したgrading HはQ-skewで `[H,M]=M` を満たすので、
N(M)はskew centralizerとこのgrading一次元の和である。

## 3. integral Sが許すcentralizer次元

どの固有値λについても、`q2−λq1` のrankが二以下なら、
複素数体上でlinear factorsの積又は平方となる。
prime I_Sがそのquadricを含むためlinear factorも含み、Sの非退化性に反する。
従って各λのgeometric multiplicityは二以下である。

Q-self-adjoint Jordan decompositionでは各Jordan blockをQ-orthogonalに取り、
Qを各blockのanti-diagonal formへ直せる。
この標準形は、最大Jordan長mに対してQ(v,M^(m−1)v)≠0となるvを取り、
非退化cyclic subspaceを分離してorthogonal complementへ帰納することで得られる。
cyclic block内のcommutantはMのpolynomialで全てself-adjointなので、skew部分は零。
同一固有値の二block長s,t間のintertwinerはdimension min(s,t)、
反対向きのblockはQ-skew条件で一意に決まる。

よってskew centralizerの次元cは、二blockを持つ各固有値についてmin(s,t)を足したもの。
全dimension五とgeometric multiplicity≤2から `c≤2`。
従って

`k=dim N(M)=c (r≥2), c+1 (r=1)`、全必要sourceのaffine dimensionは2k以下。

以下では、この次元評価よりも直接なcommon-image factorを使ってr≥2を全て除き、
r=1の三Jordan型だけを残す。

## 4. r≥2では全Dにfixed linear factorがある

この場合、全A/Bはskew centralizerに入る。
同じ固有値の二Jordan block長をs≥tとすると、全off-diagonal intertwinerのimageは

`ker(N_big^t)⊕small_block`

という固定二t次元空間へ入る。
一blockのskew centralizerは零で、異なる固有値間のintertwinerも零である。
従って全skew centralizerのimageは、dimension `Σ2min(s,t)=2c≤4` の固定空間Wに入る。
A=T Q^-1なのでAとTのimageは同じで、全ΠのimageもWに入る。

dimW≤3ならΠ²=0。
dimW=4ならWをannihilateする非零linear functional ellに対し、
実D=Euler∧Π²は `ell(x)*quartic` となる。
従ってr≥2の全Jordan型は、固有値間のcrossratioを含めて、既約被約quinticを持たない。

## 5. r=1の全三型と、唯一残る六sourceの最弱上界

scalarを引いたMはnilpotentで、geometric multiplicity≤2からpartitionは
5、4+1、3+2だけである。
Jordan5ではc=0、N(M)はgrading一次元だけなので、Πはquadric*fixed skew matrix。
実Dはquadric²*linearである。

Jordan4+1ではc=1、そのskew centralizer generator Rのimageは二次元なのでbivector C=R Q^-1はdecomposable。
grading SとRに対し、全Πは `alpha(q1,q2)*S Q^-1+beta(q1,q2)*R Q^-1`。
C∧C=0なのでΠ∧Πはquadratic alphaで割れ、実Dも同じquadratic factorを持つ。
matrixとしてR²=0と主張しているのではない。
従ってこの型も既約quinticではない。

唯一残るJordan3+2ではc=2、k=3、全必要linear sourceは六次元。
次の明示的な四つのlinear automorphismが、そのpencilのpreserving groupに入る。
Rを五×五のsuperdiagonalが一のnilpotent matrixとし、
Qをanti-diagonal entries(+1,−1,+1,−1,+1)のsymmetric matrixとする。
`N=R²` はself-adjointでJordan3+2となり、

`S=diag(1,1/2,0,−1/2,−1)`

はQ-skew、`[S,N]=N`、rankS=4である。
skew centralizerはR,R³をbasisとする。
さらに `K=(NS+SN)/2` と置くとKはQ-skew、[K,N]=N²である。
`L=N/2+K` というlinear fieldは直接

`delta_L q1=q2, delta_L q2=0`

を満たす。一方Sはdelta q1=0、delta q2=−q2、R,R³は両quadricを固定する。
S,R,R³,Lはprojective scalarを商しても独立である。
実際、そのlinear結合がscalarならQ-self-adjoint partからLの係数とscalarが零、
NとのcommutatorからSの係数も零となり、R,R³の独立性で残りも零となる。
従って実際のalgebraic preserving groupはdimension四以上である。
未証明のMöbius liftをこの下界に使う必要は無い。

Poisson条件が六次元linear sourceのproper閉条件であることを次に示す。
それなら全Poisson源≤5、projective像≤ `24+5−1−4=24` が従う。
A=S Q^-1はrank四のconstant skew matrixで、Π=q1 Aは上の必要linear sourceに入る。
しかしこのΠはprojective Poissonではない。

その一般的理由は次の通りである。
rank四constant Aのkernel coordinateをx4とし、残り四次元でAをsymplectic標準形に直す。
`[qA,qA]` は非零定数を除き `q*A^sharp(dq)∧A`。
Eulerとのwedgeの、三symplectic方向とx4を持つ成分が零なら
`x4*q*(A^sharp(dq)∧A)=0`。
四次元symplectic空間ではv→v∧Aは同型なので、A^sharp(dq)=0。
従ってqはkernel coordinateの平方だけとなるが、q1は非退化であり矛盾する。
よって少なくとも一つのJacobi polynomialは六次元linear source上で非零である。
全sourceはirreducible affine六空間なので、そのPoisson閉部分の全componentはdimension五以下。
一つの非Poisson点を全支持に一般化したのではなく、proper閉部分の次元上界という十分な矢印だけを使っている。

単一固有値の位置はpencilのlinear変更で零へ移せ、Jordan標準形に他の連続parameterは無い。
従ってこの3+2型のSは一つのPGL5 orbitであり、上のfixed-source boundが全型のPGL飽和に適用する。
これで全非退化member pencilに対する命題が従う。

## 6. 保持する補集合

pencilの全memberが退化する場合にはQ^-1が存在せず、本証明は使わない。
Kronecker singular blockを持つintegral非normal CI22や、common radicalを持つelliptic quartic coneを残す。
本命題は全rank四点の不存在又は各sourceのE帰属ではなく、既約quintic枝のG非優越である。
