# rank零plane・m二・rankM二・M(p)非零：全first-symbol family

2026-11-07、主義務G。producer幾何還元。higher Jacobiのdimension boundはまだ求めていない。
O(1) co-Higgs型、rankM=2、M(p)≠0の場合を全て保持する。
前掲 `PLANE_RANK_ZERO_MULTIPLICITY_TWO_FIRST_JET.md` の記法を使い、
normal二変数n=(z,w)、q_A=det(n,A n)、H=diag(1,−1)、E=[[0,1],[0,0]]、F=[[0,0],[1,0]]とする。

## 1. 全Lax式の一般解の有限幾何分岐

Uはconstant base vector、p=[U]、M:C³→sl2はrank二、M0=M(U)≠0。
全weight零Jacobiは

`T(M)−tM+[R,M]=0`, `[T,U]=−tU`, `R=S+M0/2`, `K=q_S`

であり、特に[R,M0]=0。
従って `R=τM0`。
さらにimM=Wについてad_R(W)⊂Wである。
Wの二whole orbitはnondegenerate plane <F,E>とBorel plane <H,E>。
nondegenerate WのnormalizerはCHであり、W内のnonzero M0のcentralizerと交わらないので、
この場合はR=0を強制する。
Borel WではR=τM0の全τが可能で、τ=0も捨てない。

kerMのpoint kはpと異なる。
baseのkernel-coordinate shearでpのk方向成分を消せる。
normal conjugation、tensor scalar、baseの可逆basis changeで、
W及びM0のsemisimple/nilpotentを下のstandard形へ直せる。
normal GL2とbase changeの逆を戻せば元の全sourceへ戻る。
M0が非零であるため必要なscalar以外をinvertしていない。

## 2. nondegenerate W：四whole family

M=bF+cEはaに依存しない。
M0 semisimpleならM0=F+E、U=∂b+∂c、K=−(z²−w²)/2と取れる。
M0 nilpotentならM0=F、U=∂b、K=−z²/2と取れる。
R=0なのでLaxは、T−tIのimageがkerM=<∂a>へ入ることを言う。
traceT=0とTp=t pも課すと全Tは

`M0 ss: T=2s*a∂a−s*b∂b−s*c∂c+ℓ(b−c)∂a`、
`M0 nil: T=2s*a∂a−s*b∂b−s*c∂c+ℓ*c∂a`。

s≠0ならactual a-shearでℓを消せる。
M0 ssではsは全nonzero scalar parameterとして保持する。
M0 nilではtensor scalar σ、normal conjugation F→σ^−1F、E→σE及び
base `a'=a,b'=b,c'=σ²c` によりM0とMを保存したままs→σsと出来るのでs=1へ直せる。
s=0ならm_D=2はℓ≠0を強制し、a-scalingでℓ=1へ直す。
各逆変換は非零scalarの逆とshearの符号反転である。

従って次の全四sourceとなる。

| 名 | U | K | T |
|---|---|---|---|
| NDSS(s), s≠0 | ∂b+∂c | −(z²−w²)/2 | 2sa∂a−sb∂b−sc∂c |
| NDSSN | ∂b+∂c | −(z²−w²)/2 | (b−c)∂a |
| NDNSS | ∂b | −z²/2 | 2a∂a−b∂b−c∂c |
| NDN | ∂b | −z²/2 | c∂a |

s=0,ℓ=0はleading normal Pfaffianが零でm≥3となり、本stratumの外として残す。

## 3. Borel W、M0 semisimple：二parameter openとresonant row

M=bH+cE、M0=H、U=∂bと取る。
全Lax solutionは

`K=(1−2τ)zw`、
`T=(2τ−2s)a∂a+s b∂b+(s−2τ)c∂c+ℓ c∂a`。

Δ=4τ−3sと置く。
actual shear a'=a+v cはℓ→ℓ−vΔとするので、Δ≠0でℓを消せる。
この場合のleading normal Pfaffianは非零constantを除いて

`Δ*a*c*q_M(n)`。

従ってrow零のΔ=0全体はm≥3であり、m二sourceとして採用しない。
Δ=0でもℓ≠0ならleading base factorはℓ*c²となりm二。
a-scalingでℓ=1に出来るので、これを別のwhole one-parameter familyとして保持する。

全m二sourceは厳密に

`BSS(τ,s): Δ≠0, ℓ=0`（全二parameter）、
`BSSR(τ): s=4τ/3, ℓ=1`（全一parameter）

となる。
τ=0は両方に存在し得るので除かない。
BSSR(0)もT=c∂aを持つm二sourceである。
Δで割った正規化の補集合を落とさず、row零のm≥3とrow非零のm二を区別した。

## 4. Borel W、M0 nilpotent：全二つの一parameter源

M=bH+cE、M0=E、U=∂cと取る。
全Lax solutionは

`K=(1/2−τ)w²`、
`T=−2s*a∂a+s*b∂b+(s*c+2τb)∂c+ℓb∂a`。

s≠0ならa'=a+v bでℓを消し、tensor scalar σとnormal conjugation E→σ^−1E、
base `a'=a,b'=σb,c'=c` を組み合わせることで、M0とMを保存してs→σs、τはそのままとなる。
σ=1/sを選びs=1へ直す。
s=0ならTとUの独立性はℓ≠0に等しく、a-scalingでℓ=1へ直せる。
ℓ=0ではT=2τb Uとなりleading Pfaffianが零、m≥3なので外に保持する。

従って全二源は

`BNS(τ): T=−2a∂a+b∂b+(c+2τb)∂c`、
`BNR(τ): T=b∂a+2τb∂c`

であり、τはいずれも全Cを走る。
必要なzero/nonzero complement及びnormal/base/tensorの実逆は以上で明示した。

## 5. 全affine symbol base次元とoperator目標

fixed planeでp二次元、rank二M八次元からpair baseは十次元。
WがBorelという条件で一、M0 nilpotentという条件でさらに一を引く。
R=0のnondegenerate WではTの全二parameterを加える。
Borelではτ一とTの全二parameterを加える。
resonance/row conditionは次の通り全source次元へ戻る。

| source | affine symbol dimension | 追加parameterの特別化 |
|---|---:|---|
| NDSS(s), s≠0 | 12 | fixed sなら11 |
| NDSSN | 11 | — |
| NDNSS | 11 | — |
| NDN | 10 | — |
| BSS(τ,s), Δ≠0 | 12 | divisor wallで高々11、finite点で高々10 |
| BSSR(τ) | 11 | fixedτなら10 |
| BNS(τ) | 11 | fixedτなら10 |
| BNR(τ) | 10 | fixedτなら9 |

plane choice六、tensor projectivization一を足し引きする。
例えばbase十二のgeneric sourceにはhigher-fiber d≤7、codimension一wallにはd≤8、
codimension二にはd≤9がG≤24に十分である。
全30/9 linearized operatorでこの弱いboundを得ることが次の有限課題である。

二parameter sourceでは、chosen maximal minorの**全**irreducible factor及びintersectionを保持する必要がある。
linear wall arrangementなら各generic lineを一parameterで、全intersectionをexactに検査出来る。
nonlinear wallが現れればexact equationを保存して停止し、未検査の曲線をgenericとして捨てない。
このnoteはまだkernelを計算せず、全G又はこのstratumの非優越を証明したとはしない。
