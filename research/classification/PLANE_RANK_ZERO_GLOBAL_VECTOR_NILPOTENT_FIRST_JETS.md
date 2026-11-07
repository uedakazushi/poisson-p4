# rank零plane・m二・globalV・normal冪零型の全first-symbol族

2026-11-07、主義務G。producerによる全sourceの手還元。
本noteはhigher-weight operatorの計算又は全Gの完了を主張しない。
既存のimmutable runner、証明書、first-jet noteを変更しない。

## 1. 記法と全weight零Jacobi

`PLANE_RANK_ZERO_MULTIPLICITY_TWO_FIRST_JET.md` §1のglobal-vector-field
supersourceを扱う。planeはP=(z=w=0)、π|P=0、実際の被約既約quintic Dの
Pに沿う一般位数は二とする。Vのdivisorial zeroを許すので、O(1)係数mapのrank一も含む。
constant normal sl2 matrixが非零冪零の場合、全weight零sourceはnormal GL2で

`Π0=T∧(z∂z+w∂w)+V∧z∂w+(a z²+bzw+c w²)∂z∧∂w`

と書ける。T,Vはtrace-freeな三変数linear vector fieldである。
base関数f,z,wに対する全Jacobiは

`z²([T,V]+bV)f+2c zw Vf`

となる。したがってV≠0からc=0、`[T,V]=−bV`。
base二関数を含むJacobiは自動であり、この二条件がweight零の必要十分条件である。
projective Euler wedgeにより隠れるscalar matrixはtrace零から零となる。
実際のleading normal quadratic Pfaffianが非零である条件は

`det(Euler_base,T,V)≠0`。                                      (1)

特にT,Vは独立である。a,bのいずれも最初から固定しない。

以下のEijは**係数行列**のmatrix unitを表し、対応するfieldはxiではなく
`xj∂xi`である。従ってfield bracketは係数行列のcommutatorの負。
actual base change gのpushforwardは係数行列を`A→gAg^−1`へ送る。
この符号の区別を全normal formで固定する。

## 2. b非零の全normalization

normal shear `w'=w+(a/b)z` はEuler_normalとz∂wを保ち、Kのaを零へ送る。
逆はminusのshearである。tensor scalar 1/bによりb=1に直し、

`K=zw`, `[T,V]=−V`

を得る。各操作はactual可逆linear change又は非零tensor scalarであり、
全higher係数も同時に変換する。chosen first symbol以外を捨てない。

field bracket条件を係数行列で書くと`TV−VT=V`である。
ゆえに`[T,V^j]_matrix=jV^j`、traceを取ってtr(V^j)=0（j=1,2,3）。
Vの全固有値は零なので、V≠0はJordan3又はJordan21だけで尽くされる。

### 2.1. Jordan3は一つのwhole orbit

`J=E01+E12`, `G=diag(1,0,−1)`とする。
Vはactual GL3 conjugationでJへ直る。
`[G,J]_matrix=J`なので全Tは

`T=G+xJ+yJ²`

である。centralizer内のactual change `g=exp(xJ+(y/2)J²)`は
V=Jを保ちTをGへ送る。指数は冪零なので有限polynomialで、逆は指数の符号反転。
したがって全sourceの代表は

`T=G`, `V=J`, `K=zw`。                                      (R)

`det(Euler_base,G,J)=x2(2x0x2−x1²)`は非零であり、(1)を満たす。

### 2.2. Jordan21は一parameter族と二つの追加orbit

`N=E01`, `B=E02`, `C=E21`, `H=diag(1,1,−2)`とする。
V=Nに直すと、全Tは

`T=Dρ+xN+yB+zC`, `Dρ=diag(ρ+1,ρ,−1−2ρ)`。                  (2)

これは全九係数へ`[T,N]_matrix=N`, trT=0を課した解で、欠けたrowはない。
Nを保つcentralizer内のunipotent changesを使う。
`Ad(exp(sB))T`ではyが`y−s(2+3ρ)`へ変わり、zは変わらず、xだけが追加で変わる。
`Ad(exp(tC))T`ではzが`z+t(1+3ρ)`へ変わり、yは変わらず、xだけが追加で変わる。
最後に`Ad(exp(qN))T`はxをx−qへ送る。
全逆はs,t,qの符号を変えた有限unipotent changeである。

したがってρ≠−2/3,−1/3ではy,z,xをすべて消去出来る。
ρ=−2/3ではz,xを消した後yだけが残り、
ρ=−1/3ではy,xを消した後zだけが残る。
Nを保つdiagonal change diag(r,r,t)はBを(r/t)B、Cを(t/r)Cへ送るので、
それぞれの残った非零係数はactualに1へ直る。零の場合も保持する。

全(2)について

`det(Euler_base,T,N)=x1²*(z*x1−(1+3ρ)*x2)`。                (3)

ゆえにρ=−1/3,z=0はm二ではなくm三以上である。
この補集合を明記して、m二の全sourceは次の三族で尽くされる。

| model | T | V | K | 全parameter |
|---|---|---|---|---|
| NDIAG(ρ) | Dρ | N | zw | ρ∈C、ρ≠−1/3 |
| NB | D(−2/3)+B | N | zw | 一つのwhole orbit |
| NC | D(−1/3)+C | N | zw | 一つのwhole orbit |

NDIAGはρ=−2/3のy=0も含む。NBとNCはresonant非零rowの全orbitであり、
標本係数の選択ではない。式(3)によりいずれもm二である。

## 3. b零ではVの全Jordan型を残す

b=0の場合、全条件は`[T,V]=0`、`K=a z²`である。
**Vを冪零とは限らない。** b非零のtrace argumentはここへ適用出来ない。
`PLANE_RANK_ZERO_GLOBAL_VECTOR_SEMISIMPLE_FIRST_JETS.md` §2の全abelian
two-plane分類を適用し、A=span(T,V)はCartan、mixed、regular-nilcentralizer、
common-kernelの四型をすべて保持する。
common-image型だけはpointwise proportionalにより(1)から除外される。

この分岐ではactual normal diagonal changeがz∂wとaを同じ非零scalar λでscaleする。
tensor scalar σも使うと

`(T,V,a)→(σT,ηV,ηa)`, `η=σλ`

となり、σ,ηは独立な非零scalarを走る。
以下のnormal formsはこの独立な二scalarとactual base normalizerで得られる。
aをtensorから独立に任意固定したものではない。

### 3.1. Cartanの全三parameter族

simultaneously diagonalな独立T,Vについて、両方のeigenvalue差が非零となる同じ座標対が存在する。
三座標対のうちnonconstantな各diagonal fieldで差が零となるのは高々一対なので、
二つの零対の和は三対を覆わないからである。
permutationと独立scalars σ,ηを使い、その差をともに1へ直す。
scalar Eulerを差し引いてtrace-free liftへ戻すと、全sourceは

`T=diag(0,1,t)_tf`, `V=diag(0,1,s)_tf`, `K=k z²`, `s≠t`

で覆われる。t,s,kは全Cを走り、有限permutationの逆と非零差の逆で元へ戻る。
leading D0は非零scalar倍の`x0*x1*x2*z²`であり、SNCではない。
normal半単純Cartanの場合のSNC component認識をこのsourceへ流用しない。

### 3.2. mixedの全有限parameter族

`H=diag(1,1,−2)`, `N=E01`とする。
actual base normalizerはH固定、N任意非零scaleを実現する。
`T=A H+B N`, `V=C H+D N`, `AD−BC≠0`の全係数を保持して、
Vのsemisimple部分とnilpotent部分の零非零を分ける。
V=H+N、V=H、V=Nの三場合へ直し、その後の残る独立T-scalarで正規化すると次の全族を得る。

| T | V | K | 全parameter |
|---|---|---|---|
| H+rN | H+N | kz² | r≠1、k∈C |
| N | H+N | kz² | k∈C |
| H+N | H | kz² | k∈C |
| N | H | kz² | k∈C |
| H+N | N | kz² | k∈C |
| H | N | εz² | ε=0又は1 |

最初の二行ではVの二係数が非零、次の二行ではVがsemisimple、最後の二行ではVがnilpotent。
例えば最後の二行の前段T=H+uN,V=Nでは、N-scale λとη=1/λがVを戻し、
uをλu、aをa/λへ送る。u≠0ならλ=1/uでT=H+Nとなり、aは全parameterのまま残る。
u=0ならaの零非零を分けてεへ直せる。
他の行は非零H係数又はN係数をσ,η及びN-scaleで割るだけである。
各逆変換は同じ非零scalarsの逆であり、零係数の場合を表の別行に保持する。

### 3.3. regular-nilcentralizerの全族

`J=E01+E12`とする。actual base normalizerは

`J→λJ+μJ²`, `J²→λ²J²`, `λ≠0`

を実現する。inverse polynomial substitutionはλ^−1、−μ/λ³。
VがregularならV=Jへ直しT=uJ+J²と書ける。
Vを固定するλとη=1/λ、及びTのJ²係数を戻すσ=1/λ²により
`u→u/λ`, `a→a/λ`となる。
u非零と零を分け、Vがsquareの場合も残すと全sourceは

| T | V | K | 全parameter |
|---|---|---|---|
| J+J² | J | kz² | k∈C |
| J² | J | εz² | ε=0又は1 |
| J | J² | εz² | ε=0又は1 |

となる。最後の行では残るλがaをa/λ²へ送るので、C上で全非零aを1へ直せる。
square rootのいずれのchoiceも可逆actual changeであり、a=0は別orbitに保持する。

### 3.4. common-kernelの二orbit

`A=<E02,E12>`のordered independent pairをactual image-plane GL2で
`T=E02`, `V=E12`へ直す。
base diagonalと独立σ,ηを組み合わせてT,Vを固定したままaを任意非零scale出来る。
従って全sourceは`K=0`又は`K=z²`の二whole orbitである。
`det(Euler_base,T,V)=x2³`は非零なので、両者ともm二である。

## 4. 全first-symbol base dimensionと必要な次の証明

b非零では、normal nilpotent **line**のchoiceが一次元である。
Vのnonzero scalarはVそのもののorbitに既に含め、二重に数えない。
regular-nil VのPGL3 orbitは六次元、Tの解空間はcentralizer二次元、Kのa,bは二次元。
rank一nil Vではそれぞれ四次元、centralizer四次元となる。
従って両whole sourcesのfixed-plane affine symbol dimensionは十一。
NDIAGのρ固定、及びNB/NCは少なくとも一つのinvariant条件を課すので上界十である。
ρ=−2/3でresonant row零という小さい層をこの上界十からさらに縮める必要はない。

b零ではAのPGL3 orbit dimensionはCartan六、mixed五、regular-nil四、common-kernel二。
ordered pair(T,V)四、normal nilpotent line一、Kのa一を加えると

| b=0 whole A source | fixed-plane affine symbol dimension upper |
|---|---:|
| Cartan | 12 |
| mixed | 11 |
| regular-nilcentralizer | 10 |
| common-kernel | 8 |

である。a零又はpairのspecial Jordan層はこのsupersource上界の中に保持される。
全higher weight一30、weight二9に対するkernel upperの和をdとすれば、
plane choice六とprojectivization一から全像dimensionは各表の値にd+5を加えたもの。
従ってb非零全sourceではd≤8、ρ固定又は追加orbitではd≤9がG≤24に十分。
b零のCartan/mixed/regular-nil/common-kernelでは順にd≤7/8/9/11が十分である。
operator minorのparameter wallではその実base dimensionを戻す必要がある。

本noteはactual全normalizationと有限parameter族、正しい補集合、dimension targetまでを証明する。
全30/9 Jacobi operatorのrank及び全wallをまだ証明していない。
m_D≥3、scalar co-Higgs first symbol、normal半単純の非Cartan higher sourceは別の未解決義務。
normal nilpotentのb=0 Cartanにも、SNC endpointを実証せずにMPS帰属を主張しない。
