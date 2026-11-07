# rank零plane・multiplicity二の全first jetとfull-rank係数map

2026-11-07、主義務G。以下は全source還元であり、全plane枝のG完了とはしない。
P=(z=w=0)⊂P4、π|P=0、実Dは既約被約quintic、Pに沿う一般位数が二とする。
位数三・四及び非Frobenius first jetは明示補集合として保持する。

## 1. 全normal weightとco-Higgs線subsheaf

P上の零性を満たす全projective sourceはdimension80（I=(z,w)の全90二次raw係数、Euler十冗長性）。
normal weightをz,wに一、normal derivativeに−一とする。
Frobenius first-jet補題からconormal idealはabelian、従ってnormal-normalのlinear項六係数は零。
残る全sourceはweight零35、weight一30、weight二9である。

weight零のbase-normal第一symbolは
`θ∈H⁰(P,T_P)⊗End(C²)`。
Jacobiをbase二関数・normal一関数へ適用すると、θのtraceless三entryはpairwise wedge零となる。
それらの共通imageのsaturationをO(k)⊂T_Pとする。
非零global sectionがあるのでk≥0、Euler列でH⁰(T_P(−2))=0なのでk≤1。
従って全first jetは次の二種類で覆われる。

- k=0：一つのglobal linear projective field Vにconstant sl2 matrixを掛けたtraceless部分。
- k=1：pencil of linesのO(1) inclusion Uを、baseにlinearなsl2 matrix M(x)に掛けたtraceless部分。

trace部分Tは初めは全H⁰(T_P)を保持する。
k=1のinclusionはH⁰(T_P(−1))=C³の非零constant vector U、すなわちpoint p∈Pで指定される。
normal weight零のvertical quadratic K(z,w)も全三係数を保持する。
Euler gaugeのliftの選び方はKへ吸収出来る。

## 2. k=1の全weight零Jacobiは一つのLax式

base coordinatesをx=(x0,x1,x2)としU=∂x0を選ぶ。
M(x)=x0M0+x1M1+x2M2は全linear sl2 matrix、n=(z,w)、E_n=z∂z+w∂wとする。
全weight零sourceは

`Π0=T∧E_n+U∧(M(x)n·∂n)+K(n)∂z∧∂w`。           (1)

Tはtrace-free linear matrixでliftする。
`q_A(n)=det(n,A n)` と置く。全binary quadraticは一意にK=q_S、S∈sl2と書ける。
base二関数に関するJacobiはscalar matrix IとMの可換性から自動である。
残るJacobiのbase-vector係数は

`q_M[T,U]+U*(V_M K+q_T(M)+det(Mn,M0n))`、           (2)

ここでV_M=Mn·∂n。
projective Jacobiは(2)とEuler_baseのwedgeが零という正確な条件。
もし[T,U]がUと独立なら、(2)は固定二次元base-vector spaceに入り、
全xでEuler_baseにproportionalであることは不可能なので(2)=0となる。
その二成分を比べるとq_M=0、M=0となり非零traceless部分に反する。
よってTはpを保存し、[T,U]=−t00 Uである。

二×二traceless行列の恒等式

`V_M q_S=q_[S,M]`, `det(Mn,M0n)=q_[M0,M]/2`

により、(2)は厳密に

`T(M)−t00 M+[R,M]=0`, `R=S+M0/2`。                (3)

x0係数は[R,M0]=0である。
任意のrankのMについて(1)–(3)は全weight零必要十分条件であり、chosen pointだけの診断ではない。

## 3. rankM=3の全実normalizationと二つの一parameter源

ここからM:C³→sl2がisomorphismの場合に限定し、他のrankを保持する。
(3)をmatrix式にすると

`T=t00 I−M^−1 ad_R M`。

traceT=0からt00=0。
M0=M(p)は非零、centralizer_sl2(M0)は一次元なのでR=τM0となる。
従って

`T=−τ M^−1 ad_M0 M`, `K=q_(τ−1/2)M0`。          (4)

τ=0ならT=0でleading normal quadratic Pfaffianが零となり、m_D=2に反する。
従ってこのstratumではτ≠0。
baseの可逆linear transformationでMをstandard sl2 coordinatesへ直す。
normal GL2 conjugationとtensor scalarでM(p)をsemisimple又はnonzero nilpotentに直せる。
逆変換はMのinverse、normal GL2 inverse及び非零scalarであり、元の全sourceへ戻せる。
normalizationでunknown higher-weight係数を制約しない。

standard M(a,b,c)=[[b,c],[a,−b]]に対し全二源は

| 型 | U | T | K |
|---|---|---|---|
| semisimple | ∂b | 2τ(a∂a−c∂c) | (1−2τ)zw |
| nilpotent | ∂a | τ(c∂b−2b∂a) | (τ−1/2)z² |

(1)へこの全係数を代入する。
τ以外のhidden modulus又はeigenvalue complementは無い。
M0=0はrankM=3とp≠0に反する。
τ=0、rankM≤2、k=0、m_D≥3を本二源へ無理に含めない。

## 4. 実parameter次元と有限higher-weight課題

fixed planeの全rankM=3 weight零sourceは、p二次元、M九次元、τ一からaffine dimension十二。
semisimple M(p)はopenで十二、nilpotent M(p)は一つのnonzero discriminant equationで十一である。
τを一つの非零値へ固定すると各々一ずつ下がる。
plane choiceはGr(3,5)のdimension六、tensor projectivizationは一を引く。
fixed normalized Π0のhigher-weight full solution fiberをdとすると、全PGL像上界は

`17+d`（semisimple一般τ）、`16+d`（nilpotent一般τ）。

finite τ-wallに固定した層ではさらに一を引く。
従ってgenericでd≤7/8、finite wallでd≤8/9がG≤24に十分である。
これはfixed-symbol sliceのdimensionをそのままambient像へ流用したものではない。

全higher sourceは

`Π1=zL+wN+z²D∧∂w+zw(B∧∂z+E∧∂w)+w²C∧∂z`（全30）、
`Π2=Σ_basepairs(αz²+βzw+γw²)∂i∧∂j`（全9）

である。
weight一Jacobiは全[Π0,Π1]、fixed Π1に対するweight二のsolution fiberは
全[Π0,Π2]のkernelのtranslateなので、
`d≤dimkerD1+dimkerD2` である。
全350 projective coefficientを使い、一parameterのnonzero minorの全例外根を保持する有限jobを行う。
generic kernelが実shear六だけでD2 injectiveなら、
base change x→x+C(z,w)はactual finite unipotent linear transformationでΠ1を消し、Π2も消える。
そのときDはΠ0のbase-quadratic因子を持つので既約quintic枝は空。
special τでは上の弱いdimension boundだけで十分である。

この有限jobは30/9 linear operatorであり、全80-source nonlinear GBを求めるものではない。
全plane G又は未処理の他first-jet strataの完了は主張しない。
