# rank零plane・m二・rank M二・M(p)=0の独立review

2026-11-07。主義務G。producerと別作者による全source幾何の読解及び独立Fraction算術。**PASS**。

対象はrank零plane Pを持つPoisson tensorで、実既約被約quintic DのPに沿う一般位数が二、
coHiggs共通imageのsaturationがO(1)、linear coefficient map M:C³→sl₂のrankが二、
そのkernel pointがO(1) inclusionのpoint pと一致する場合である。
**この全stratumのprojective像dimensionは二十二以下**。
既知成分への全点帰属、非線形全solution、他のplane枝又は全G完了を主張しない。

## 1. 幾何入口と全四正常形

`PLANE_RANK_ZERO_MULTIPLICITY_TWO_FIRST_JET.md` のm二/Frobenius入口を読解した。
rank二共変idealはabelian、traceless三vector fieldのpairwise wedgeは零となる。
共通imageをsaturateしたrank一subsheafは、smooth P²上でtorsionfree quotientを持つためreflexive、
従ってO(k)である。非零global sectionとEuler列のH⁰(TP²(−2))=0からk=0又は1となる。
以下はk=1だけを扱う。

U=∂a、M=bM₁+cM₂とし、normal variablesをn=(z,w)、Eₙ=z∂z+w∂wとする。
全weight零sourceは

`Π₀=T∧Eₙ+U∧(Mn·∂n)+q_R(n)∂z∧∂w`, `q_R=det(n,Rn)`

である。ここではM(p)=0なので一般式のM₀/2補正は零である。
tracefree Tを用いた全Lax式は

`R∈N_sl₂(im M)`, `A=tI−J`, `3t=tr J`,

`T=(ta+r_b b+r_c c)∂a+A(b,c)·∂(b,c)`, `J=ad_R|im M`。
base shear a'=a+ℓ(b,c)はrow rをr−ℓJへ変える。
R=0ならΠ₀はUを共通因子に持ちΠ₀²=0となる。これはDのleading normal二次symbolを零にし、
m二に反する。したがってR=0を黙って除外したものではない。

Killing orthogonal lineによりsl₂内の二planeは二whole orbitだけである。
H=diag(1,−1)、E=E₁₂、F=E₂₁とすると、
nonisotropic型W=〈F,E〉のnormalizerはCH、isotropic型W=〈H,E〉のnormalizerは〈H,E〉。
これは[H,E]=2E、[H,F]=−2F、[E,F]=Hから直接確認できる。
M₁,M₂は独立なので、各Wのbasisをbase changeで直せる。

この全sourceは次の四正常形へ戻る。

|型|M|T|K=q_R|
|---|---|---|---|
|ND|bF+cE|2b∂b−2c∂c|−2zw|
|BSSρ|bH+cE|(2a+ρb)∂a+2b∂b−4c∂c|−6zw|
|BN1|bH+cE|c∂a+2b∂c|−w²|
|BN0|bH+cE|2b∂c|−w²|

BSSでは**全ρ∈C**を保つ。nilpotent rowだけを零・非零へ分ける。

## 2. 実可逆normalizationの全parameter恒等式

独立checkerは、sample値ではなくα,β,r_b,r_c,εを残したLaurent環で、
以下の全変換、両方向のinverse composition、全十bivector coefficientの一致を確認した。
全変換はbase/normal gradingを保存し、higher coefficientへの制限を導入しない。

NDではR=αH、α≠0、T=(r_b b+r_c c)∂a+2αb∂b−2αc∂cに対し

`a'=a−r_b b/(2α)+r_c c/(2α), b'=b/α, c'=c/α, z'=z, w'=w`,

tensor scalar1/αでNDになる。逆は
`a=a'+r_b b'/2−r_c c'/2, b=αb', c=αc'`。

BSSではR=αH+βE、α≠0、
T=(2αa/3+r_b b+r_c c)∂a+(2αb/3)∂b+(2βb−4αc/3)∂cに対し

`a'=a+r_c c/(2α)−r_cβb/(2α²), b'=3b/α, c'=3c/α−3βb/α²`,

`z'=z+βw/(2α), w'=w`, tensor scalar3/α。
残るparameterは`ρ=r_b+βr_c/α`であり、逆は

`a=a'−r_c c'/6, b=αb'/3, c=αc'/3+βb'/3, z=z'−βw'/(2α), w=w'`。

Borelのα=0、β≠0はnilpotent chartとして保持する。
R=βE、T=(r_b b+r_c c)∂a+2βb∂cに対し

`a'=a−r_b c/(2β), b'=b, c'=c/β, z'=z/β, w'=w`

でBNε、ε=βr_cとなる。逆は
`a=a'+r_b c'/2, c=βc', z=βz'`、b,wはそのままである。
ε=0はBN0。ε≠0では

`a'=a, b'=b/ε, c'=c, z'=εz, w'=w`, tensor scalar1/ε

がBN1へ戻し、逆はb=εb',z=z'/εである。
GL₂のscalarを許すので、この書き方には平方根選択も不要である。

## 3. 全30/9 sourceと全350行

独立runnerはproducerをimportせず、既存の独立Fraction sparse arithmeticと
標準ライブラリによるQQ(ρ)演算だけを再利用する。
全weight一raw方向三十六に対し、Euler gauge六方向とcanonical三十方向の合併rankが三十六である。
その六gaugeはE∧(zY+wZ)、Y,Z∈C³であり、全weight一sourceを尽くす。
weight二はraw九方向を全て保ち、対応するlinear Euler gaugeは存在しない。
reportに保存された全independent symbolic parameterもこの順序で直接抽出して照合した。

各Π₀からprojective Schoutenの全350係数を独立展開し、
350×30のD₁、350×9のD₂を全entryで照合した。
actual linear shear x→x+C(z,w)のpushforward微分を構成し、Euler補正でcanonical三十方向へ戻して
30×6行列も全entryで照合した。chosen slice又は一部のJacobi式への置換はない。

QQ(ρ)上ではcomplete kernelの全行代入、kernel独立性、nonzero maximal minorを確認した。
D₁,D₂及びshearのentryが全てρの多項式であることも確認した。
BSS selected minor productは

`6112251732438835879213006848 (ρ−3)(ρ−1)²(ρ+1)²(ρ+3)`。

独立univariate除算で、この四root以外に因子がないことを確認した。
**ρ=−3,−1,1,3の全て**で全行列をexact QQへspecializeし、
complete kernel及び新しいnonzero maximal minorを再確認した。ρ=0はgeneric minorが非零である。
したがってgeneric rankだけから境界を捨てていない。

|型|dim ker D₁|dim ker D₂|actual shear rank|
|---|---:|---:|---:|
|ND|6|0|6|
|BSS genericρ|6|0|6|
|BSS ρ=3,1,−1|6|0|6|
|BSS ρ=−3|7|0|6|
|BN1|6|0|6|
|BN0|7|0|6|

Π₀自身の全projective Jacobiと実divisorも直接展開した。
D₀はNDで4bc(bz²−cw²)、BSSで−6bcw(2bz+cw)、BNで2b²w(2bz+cw)となる。
これはsymbolの実divisorであり、full higher sourceの因子型をこれだけから結論しない。

## 4. whole sourceへのdimension復元

fixed planeでpのchoiceは二、pをannihilateするrank二map Mのchoiceは六、合計八である。
NDはopenで、normalizer R一及びrow二を足して十一。
Borel image-plane条件は一式であるためM部分は七となる。
BSSではR二とrow二を足して十一、BNではnonzero R一とrow二を足して十。
BN0のε=0はrowの一つの線形条件なので九である。

BSSのfixedρ層は十以下である。p,Mを直すframeを、rowから独立に有限個のalgebraic chart上で選ぶ。
上の式ρ=r_b+βr_c/αはr_b係数一を持つから、各fixedρはrow二の中の真正一条件となる。
frameのfinite coverによるdimension増加はない。これでwall層のdimension低下を明示できる。

Π=Π₀+Π₁+Π₂の全Poisson条件のうち、weight一はD₁Π₁=0。
各固定Π₁に対し、weight二はD₂Π₂に既知のΠ₁二次式を加えたaffine linear systemである。
そのfiberは空又はker D₂のtranslateであり、高weight条件はさらに狭めるだけである。
よってhigher full fiber dimensionは`d≤dim ker D₁+dim ker D₂`。

plane choice六を足しtensor projectivization一を引くと、全像は

|whole stratum|projective image upper|
|---|---:|
|ND|11+6+6−1=22|
|BSS genericρ|11+6+6−1=22|
|BSS ρ=3,1,−1|10+6+6−1=21|
|BSS ρ=−3|10+7+6−1=22|
|BN1|10+6+6−1=21|
|BN0|9+7+6−1=21|

となる。未知dimension二十五以上の成分の一般元はこの全stratumに入らない。
M(p)≠0、rank M≤1、k=0、m_D≥3及び他のG義務は残す。

## 5. 固定証明書

有限job gate PASS後に実行し、2.282秒、最大RSS22,172KiBで終了した。
producer及びCASをimportしていない。旧rank三のproducerはASTだけを照合し、旧計算は再実行していない。
producer runnerの後続append変更には依存せず、immutable public report中のexact bytesを参照する。

|artifact|SHA256|
|---|---|
|code/classification_plane_rank_zero_rank_two_review.py|765f92b493d48c0b9d18f25d2971e6d1686d0abb4d3b8522cae3fc8d7f65dcd2|
|reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_REVIEW_20261107.json|d766209f6ae16bd90e17c878e24ad5fb3cb5c621a2ad1f74bddbf09bf5c62ef4|
|reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_REVIEW_JOB_20261107.json|af43f6afb2bb254eac7e6591315c1aabb14326700b465a7fdfbdd356b2ba33af|
|producer CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_KERNEL_WEIGHTS_20261107.json|ce13a6296477adeea413f9004fcd28e579ed4001fc40019efbb3d3a09a309d3f|
|producer PLANE_RANK_ZERO_RANK_TWO_KERNEL_FIRST_JET.md|c4de97e904a4fb1fccc3bda8164c2ce4780519e57a72fdc070d9d2ca1bc0a108|
|producer PLANE_RANK_ZERO_RANK_TWO_KERNEL_WEIGHT_BOUND.md|414336a9eb60277ac62c62f6f9b8629004a734dd4c5ee10b4798ae4329b79045|

その他の独立arithmetic helperと旧AST入力はjob及びreview reportのinput manifestへ全hashを保存した。
