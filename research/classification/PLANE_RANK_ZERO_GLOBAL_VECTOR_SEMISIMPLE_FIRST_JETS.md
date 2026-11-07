# rank零plane・m二・globalV・normal半単純型の全first-symbol族

2026-11-07、主義務G。producer幾何還元。
co-Higgs traceless部分が一つのglobal projective vector field Vとconstant nonzero semisimple sl2 matrixの積である全sourceを扱う。
Vのdivisorial zeroを許し、O(1)の係数rank一もこのsupersourceに保持する。
normal nilpotent型及びm_D≥3は別補集合である。

## 1. 全weight零Jacobi

normal basisを直すと

`Π0=z(T+V)∧∂z+w(T−V)∧∂w+kzw∂z∧∂w`, `[T,V]=0`

となる。T,Vはtrace-free linear matrixでliftする。
全K三係数からkzwだけへ制限されるJacobiは
`PLANE_RANK_ZERO_CARTAN_SNC_RECOGNITION.md` §1の全係数計算による。
m_D=2は `det(Euler_base,T,V)` が非零であることに等しい。
従ってT,Vはsl3内で独立で、その二次元span Aはabelianである。

## 2. 二次元abelian subalgebraの全四型

Aがregular matrix（minimal polynomial次数三）を含むなら、
そのsl3 centralizerは二次元であるからAそのものとなる。
regular matrixのJordan型から全三型はCartan、
mixed centralizer `<H,N>`、regular nilpotent centralizer `<J,J²>` である。
ここで

`H=diag(1,1,−2)`, `N=E01`, `J=E01+E12`。

Aがregular matrixを含まない場合を残す。
nonzero semisimple matrixがあれば、そのdouble-eigenvalue二block内のnonscalarな別のA要素と
generic linear combinationを取ることでregular matrixを得てしまう。
よって全nonzero要素はrank一nilpotentである。
rank一matrix u⊗ℓとv⊗mの和もrank一なら、u,vがdependent又はℓ,mがdependent。
従って二次元linear space全体はcommon image又はcommon kernelを持つ。
trace零からcommon-image型は `<E01,E02>`、common-kernel型は `<E02,E12>` へ直る。
前者のvector fieldsは同じconstant方向でpointwise proportionalとなり、m_D=2に反する。

従って全m二AはCartan、mixed、regular nilcentralizer、common-kernelの四型で尽くされる。
Cartanは全higher sourceがMPS40へ入ると前掲SNC noteで証明したので、以下では残る三型だけをnormal formにする。

## 3. mixed centralizerの全三族

T=aH+bN、V=cH+dN、ad−bc≠0と書く。
actual base normalizerはHを固定しNを任意の非零scalarでscaleする。
これはHのdouble-eigenvalue block内のdiagonal changeで実現し、逆も可逆である。
tensor scalar σも同時に使用する。

c,dとも非零ならσ=1/c、N-scale=c/dでV=H+Nとなり、
T=uH+vN、u−v≠0を得る。
c≠0,d=0ならV=H、T=uH+N。
c=0,d≠0ならpair independenceでa≠0、V=N、T=H+uNとなる。
各choiceは上の非零scalarの逆を戻せば元の全T,Vへ戻る。
kもtensor scalarに従い変換されるが、全Cを走るので固定して捨てない。

| model | T | V | K |
|---|---|---|---|
| MIX(u,v,k), u−v≠0 | uH+vN | H+N | kzw |
| MIXH(u,k) | uH+N | H | kzw |
| MIXN(u,k) | H+uN | N | kzw |

全三/二parameter族であって、u,v,kを標本値へ置き換えていない。

## 4. regular nilcentralizerの全二族

T=aJ+bJ²、V=cJ+dJ²、ad−bc≠0。
actual base normalizerは

`J→λJ+μJ²`, `J²→λ²J²`, λ≠0

を実現する。
右辺のregular nilpotent matrixはJと同じJordan型なので、Jordan chain basisの変更で
実際のGL3 conjugationが得られる。
逆polynomial substitutionはλ^−1及び−μ/λ³であり、degree二までで完結する。

c≠0ならtensor σ、λ=1/(σc)、μ=−dλ²/cでV=Jとなる。
TのJ係数はu=a/c、J²係数は `(bc−ad)/(σc³)`。
σ=(bc−ad)/c³を選べばT=uJ+J²となる。
c=0ならd,a≠0、λ=a/d、σ=d/a²、μ=−bλ²/aでT=J、V=J²となる。
全inverseは上の非零scalarsを戻すactual変換である。

| model | T | V | K |
|---|---|---|---|
| NILREG(u,k) | uJ+J² | J | kzw |
| NILSQUARE(k) | J | J² | kzw |

## 5. common-kernel型の全二orbit

A=<E02,E12>はcommon kernel及びimage plane <e0,e1>を持つ。
baseのこのplane上のGL2 actionはA上の任意のGL2 basis changeを実現する。
したがって独立なordered pairをT=E02、V=E12へ直せる。
tensor scalarを、base image-plane scalarの逆でT,Vについて戻すと、kだけを任意の非零scalarでscale出来る。
よってk=0、k=1の二whole orbitとなる。
逆はいずれも非零scalarの逆であり、k零の補集合を落としていない。

## 6. 全first-symbol base次元

sl3内のsubalgebra normalizerは、mixedでdimension三、regular nilcentralizerで四、
common-kernel型で六である。
mixedではHのblock centralizer内のN-line normalizer、nil型ではcentralizer二にgradingと
`[E01,J]=J²` の方向を加えることで直接求まる。
common-kernel型はそのkernel planeを保存するparabolicである。
PGL3 orbit dimensionsはそれぞれ五、四、二。

ordered basis(T,V)四、normal semisimple sl2 lineのactual orbit二、K一を加えると、
fixed-plane affine first-symbol dimensionsはmixed十二、nil十一、common-kernel九。
closed V-condition及びk零を戻して

| source | affine symbol dimension |
|---|---:|
| MIX(u,v,k), u−v≠0 | 12 |
| MIXH(u,k), MIXN(u,k) | 11 |
| NILREG(u,k) | 11 |
| NILSQUARE(k) | 10 |
| common-kernel k=1 | 9 |
| common-kernel k=0 | 8 |

を得る。
全higher fiberのlinear operator kernel upperをdとすれば、plane choice六とprojectivization一から
image upperはこの表の値にd+5を加えたもの。
例えばgeneric d=6なら順に23、22、22、21、20、19となる。
parameter wallでは実base dimensionの減少を戻し、全exceptional curve/intersectionを保持する必要がある。

全weight一30・weight二9を維持する有限jobが次の課題である。
本noteはそのrankをまだ計算しておらず、全sourceのG非優越又は全MPS帰属を未証明で主張しない。
normal nilpotent型、plane m三/四、rank二restrictionも明示残余である。
