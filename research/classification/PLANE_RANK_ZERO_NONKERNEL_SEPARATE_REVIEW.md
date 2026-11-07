# rank零plane・m二・rank M二・M(p)非零の全八family独立review

2026-11-07。主義務G。producerと別作者による全source幾何及び独立exact算術。**全八family PASS**。

対象はrank零plane P、実既約被約quintic DのPに沿う一般位数二、
共通coHiggs image O(1)、rank M=2、M(p)≠0という全stratumである。
以下は全projective像dimension≤23という弱いG結論であり、全Poisson support又はE帰属を求めない。

## 1. 全Lax式からのwhole-source分岐

kernel側の独立reviewで既に読解した全入口を用い、UをO(1) inclusion、M₀=M(U)とする。
Π₀=T∧Eₙ+U∧(Mn·∂n)+q_(R−M₀/2)∂z∧∂wについて、全weight零Jacobiは

`[T,U]=−tU`, `T(M)−tM+[R,M]=0`。

M₀≠0であり[R,M₀]=0だからR=τM₀。
W=im Mはnonisotropic orthogonal line型〈F,E〉又はisotropic型〈H,E〉の二whole orbit。
前者のnormalizer CHはW内の非零M₀のcentralizerと交わらないのでR=0となる。
後者はBorelで、τ=0も含めて全τを保つ。
ker M≠pなので、kernel-coordinate shearでpのkernel方向成分を消し、
MとM₀を以下のstandard形へ直せる。これは全可逆linear changeであり、sampleへの制限ではない。

nondegenerate WではM=bF+cE、M₀ semisimpleはU=∂b+∂c、
nilpotentはU=∂bと取る。T−tIのimageがker Mへ入ること、trace T=0、Tp=tpを合わせると、
全Tはdiag(2s,−s,−s)に、semisimpleではℓ(b−c)∂a、nilpotentではℓc∂aを加えた形である。
s≠0ならa-shearでrowを消す。semisimpleのsは全非零parameterとして保持する。
nilpotentのs≠0はnormal/base/tensorの実変換で一へ直せる。
s=0ではℓ≠0をa-scalingで一へ直す。ℓ=0ならleading normal二次divisorが零となり、m二に反する。
これがNDSS(s≠0)、NDSSN、NDNSS、NDNの全四familyである。

Borel、M₀=HではM=bH+cE、U=∂b、R=τHで、全Tは

`((2τ−2s)a+ℓc)∂a+s b∂b+(s−2τ)c∂c`。

Δ=4τ−3sについて、a'=a+v cはℓ→ℓ−vΔと作用する。
Δ≠0でrowを消したBSS(τ,s)を全二parameterで保つ。
Δ=0、ℓ=0はleading二次divisor零なのでm≥3補集合。
Δ=0、ℓ≠0はa-scalingでℓ=1へ戻るBSSR(τ)であり、全τ∈Cを別に保つ。

Borel、M₀=EではM=bH+cE、U=∂c、R=τEで、全Tは

`(−2s a+ℓb)∂a+s b∂b+(s c+2τb)∂c`。

s≠0でrowを消し、actual normal/base/tensor scaleでs=1にするBNS(τ)と、
s=0、ℓ≠0を一へ直すBNR(τ)が全sourceを尽くす。
s=ℓ=0はleading二次divisor零であり、m≥3へ送る。
τ=0を捨てず、BSSR/BNS/BNRのいずれでも直接保持する。

## 2. 実normalizationの独立全parameter検査

独立checkerは全十bivector coefficient及び両inverse compositionをLaurent環で検査した。
主な実変換は次の通りである。

- NDのrow除去はa'=a+ℓ(b−c)/(3s)又はa'=a+ℓc/(3s)。
- ND nilpotentのs-scalingは(a,b,c,z,w)→(a,b,σ²c,σz,w)、tensor scalarσ。
  M、M₀、Kをstandardのままs→σsとする。
- BSSのrow除去はa'=a+ℓc/Δ。Laurent補助変数dを用い、分母を払ってd=4τ−3sをexact代入した。
- BNSのrow除去はa'=a−ℓb/(3s)。そのs-scalingは
  (a,b,c,z,w)→(a,σb,c,z/σ,w)、tensor scalarσであり、s→σs、τは不変。
- 全四nonzero-row branchはa'=a/ℓでrow一になる。

逆はshearの符号又は非零scaleを戻すだけである。全変換はnormal gradingを保存する。
M₀ semisimpleの標準化には非零eigenvalueに依るtensor scaleを用い、そのscaleを含む元のM全parameterは§5で復元する。eigenvalueの符号choiceは有限である。
M₀ nilpotentの標準化は可逆normal conjugationであり、τを0又は1に置換するものではない。

## 3. 全30/9 operatorとexact係数体

全weight一raw三十六方向、Euler六方向、canonical三十方向、weight二の全九方向は
kernel側のcheckerを再利用してsymbolic source全parameterのまま照合した。
新算術はproducerもCASもimportしない。
全八Π₀をU,T,M,Kから別に構成し、全projective Jacobiと実divisorを再展開する。
全350×30、350×9、actual30×6 shearは全entryを保存値と比較する。
各complete kernelを全行へ代入し、free-variable identity rowsで独立性を確認し、
保存された行列式を独立に再計算する。

一parameterはFractionによるQQ(τ)。二parameter BSSはQQ[τ,s]のsparse多項式及び有理式である。
BSSのmaximal minorはfraction-free Bareiss法で計算し、各除算はremainder零を確認する。
数値specializationで二parameter恒等式を推測していない。

NDSSの全selected minor因子は、非零constantを除き

`s⁹(s−1)(s+1)(3s−1)(3s+1)(2s²−1)(3s²−1)(18s²+1)(72s²+1)`。

全rational root及び四irreducible quadratic factorを残す。
四fieldではα²=1/2,1/3,−1/18,−1/72の二次元QQ代数を独立実装し、
全matrixを代入してcomplete kernel/nonzero minorを確認する。
これらはQQ上非平方でありfieldである。全complex embeddingは忠実なので、αと−αの両rootでrankが同じである。
数値的√近似又は一方のrootだけの検査ではない。
s=0も算術上検査し、m≥3補集合として明示する。

## 4. BSS全二parameter壁と全intersection

generic QQ(τ,s) complete kernelはD₁六、D₂零。
selected minor productの全十七linear factorの積との恒等式を独立に確認する。
Δ=0 lineを除いた**全十六line**をA¹とexact同型にparametrizeし、
その全350行に対する新しいQQ(r) kernel/minorと、そのminorの全rootを検査する。
lineの部分的sampleをgeneric線全体の証明と扱っていない。

全142 line-root certificateのうち十四はΔ=0上のm≥3補集合である。
非除外の百二十八ではhigher fiber d≤8。
二parameter generic patch、各line generic patch、各line minorの全rootで全C²を覆う。
さらにcheckerは十七lineの全finite pairwise intersection百二十一組を列挙し、
各々がそのline上のnonzero minor patch、明示root certificate又はΔ=0補集合へ入ることを確認する。
intersectionをgeneric扱いして捨てていない。
Δ=0上のrow非零m二sourceは別のwhole BSSR familyとして完全に処理する。

## 5. full source次元の復元

fixed planeでp二、rank二M八を合わせて十。
W Borelは一条件、M₀ nilpotentはさらに一条件を与える。
NDは全T二parameter、Borelはτ一と全T二parameterを加える。
normalization前のsymbol dimensionは次となる。

|source|affine symbol dimension|generic d|generic projective image upper|
|---|---:|---:|---:|
|NDSS(s≠0)|12|6|23|
|NDSSN|11|6|22|
|NDNSS|11|6|22|
|NDN|10|7|22|
|BSS(τ,s),Δ≠0|12|6|23|
|BSSR(τ)|11|6|22|
|BNS(τ)|11|6|22|
|BNR(τ)|10|7|22|

canonical parameterのfinite wallはsource dimensionを一、BSS whole lineは一、finite pointは二下げる。
frameはp,M及びM₀標準化の有限algebraic chartで選べ、rowから独立に取れる。
τはR=τM₀の比であり真のscalar parameter、semisimple側のsもnormalized diagonal parameterである。
このためwallを定める非零linear式はsource上の真正codimension条件である。

fixed Π₁についてΠ₂のPoisson必要条件はD₂のinhomogeneous linear systemであり、
全higher fiber dimensionはd≤dim ker D₁+dim ker D₂。
plane choice六を加えtensor projectivization一を引く。
BSS lineは11+6+6−1=22、nonexcluded pointは10+8+6−1≤23。
一parameterの全finite wallも二十三以下であり、BNRのτ=1,3/2,1/2ではd=9だがbase九なので二十三となる。
全八familyのprojective像は従って二十三以下である。

rank M≤1、k=0、m≥3及びrank二restrictionの残余は本命題で閉じない。
Dのleading symbolの可約性をfull higher Dへ流用しない。

## 6. 固定証明書と保存範囲

finite job gate PASS後に実行し、153.124秒、最大RSS36,452KiBで全て終了した。
producer、producer helper及びCASをimportしていない。
旧kernel-side producerとの差分は全既存FunctionDef ASTを照合してappend互換性を確認しただけであり、旧計算は再実行していない。
本reviewはimmutable public report内のexact producer bytesへbindし、mutable runnerには依存しない。
全ten universal normalization及び全matrix certificateは新reportへ保存した。

|artifact|SHA256|
|---|---|
|code/classification_plane_rank_zero_nonkernel_review.py|43ce671f1a2730c48e004c6f4c33d82ed71d334f4e93b47109d6a78d58d826d3|
|reports/CLASSIFICATION_PLANE_RANK_ZERO_NONKERNEL_REVIEW_20261107.json|7006bcf9a625244f59aa10df70fec0e23d57fd8081a8bc21c6fadd83de58b386|
|reports/CLASSIFICATION_PLANE_RANK_ZERO_NONKERNEL_REVIEW_JOB_20261107.json|649c0ab410d90391ef8d5a92c08919430768c8d3fe1f64d8b3ed3d84c4d2ccdc|
|producer CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_WEIGHTS_20261107.json|c9f62505c7edfa4469296811ed9fd8f2c92f5d9ca4c6455583a28ba3737b7ce7|
|producer PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_FIRST_JET.md|027c7e8281542595300f8dcf06fde95cb7d97d686a41133678a3303dbbae6cdb|
|producer PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_WEIGHT_BOUND.md|d72e19a9af976948d185bac833f7ced8303753fc59db76b6ece60fa21eb293cb|

再利用した独立RF/sparse/Laurent arithmetic及びkernel-side checkerの全path/hashはjobとreview reportに保存した。
計数、全体台帳、原稿及びGitはintegratorが別に扱う。
