# Cubic-cone rank二の全五source：Euler枝、B全λ、A一般λの純還元

2026-11-08、PRIMARY `E_reduced_nonnormal`。既保存全27source／全五moduleと全実normalizer逆を使う純証明。
新M1、operator、rank、minor、λ生成、native、GB又はpoint評価は無い。
結論は、被約既約actualquinticの未処理sourceを **A_ss/A_nil、λ=2/5,3/5の全四fiber**へ縮めることである。
四fiberの新actualD因子／既知親／全空性は本稿では主張しない。

## 1. 全canonical sourceとactual vertex shear

`CUBIC_CONE_RANK_TWO_SEPARATE_REVIEW.md` の全源定理を入力とする。
座標(x0,x1,x2,x3,c)、E=Σ_(i<4)xi∂i、
canonical divergence-zero liftは全係数上で

`Π=P(x)+(c*d(x)+V(x))∧∂c`, `d_j=Σ_(i<4) ∂_i P_ij`。

Pは保存base九源、Vは保存全十八vertex源でdivV=0。
元350Jacobiは全base十二quadricsとhomogeneous full350×18 vertex linear systemへ正確に分かれ、
V²条件とinhomogeneousparticularを追加で省いていない。
全actualSL2／tensor scaleの逆により、全非零baseはA_ss(λ),A_nil(λ),B_nil(λ),W_ss,W_nilへ戻る。
base零はdecomposableでrank≤2。

任意のbase linear form ell に対し actual determinant-one shear `c'=c+ell(x)` は

`V -> V+Ψ(ell), Ψ(ell)=P#dell−ell*d`

を実現する。inverseは同じbase formを引くこと。
固定Sのidealはc非依存なので、whole tangency/conormal sourceを保つ。
Jacobianがconstant一なのでcanonical divergence-zero条件も保つ。
従って Ψの四列は各fullvertexkernelに入る。

## 2. A二lineの全一般fiber：四shearが全moduleを尽くす

q0=x0x2−x1²、q2=x1x3−x2²とする。
保存された全126 particular actualD coefficientsは正確に

```
A_ss: F0=5(2λ−1)*c*q0*q2,
A_nil: F0=5(2λ−1)*c*q0².
```

これらは元canonicalΠ(P,V=0)の全signedEulerPfaffianである。
λ≠1/2で Ψ(ell)=0 とすると、同じraw tensorがshearで不変となる。
一方actualDはdet1の逆substitutionにより
`5(2λ−1)(c'−ell)Q4` へ移る。
元F0との一致からell*Q4=0、Q4非零なのでell=0。
したがって四shear mapはinjective。
fullkernelが四次元の各λでは、全Vを一意なactualshearで零へ戻せる。
全sourceはproper linear factor `(c+base linear)` を持ち、対象既約quinticには入らない。

これはsourceparameter-dependentなminor chartを一つだけ選んだ逆ではない。
各λのlinear map Ψはfullcolumnrank四、target全kernelも四次元なので線形同型。
非零minorによる通常のCramer inverseの開被覆が全λを覆い、未処理のzero補集合は無い。
実shear inverse自体は全sourceでpolynomial coordinate transformationである。

## 3. λ=1/2を含む全四次元fiberのconstant-rank橋

既保存全lineと全十九wallの完全kernel theoremによると、
A_ss/A_nilのrankが十四でない有限λは、それぞれ **2/5、3/5だけ**である。
他の全minor／generic-basis denominator wallでは元350matrixへ戻してrank十四、kernel四と確認済み。
従って irreducible base U=A1\{2/5,3/5} 上のfullvertexkernelはrank四vector bundle。

λ≠1/2のopen dense部分の全sourceは§2でactualD可約又は零。
proper-factor locusはprojective multiplication mapのproper imageの有限和、
そのaffinecone（D零を含む）は閉である。
全kernelbundleのtotalspaceはirreducible、actualD coefficientmapはregular。
したがって全Uの全fiberも同じ閉条件を満たす。
**λ=1/2の全vertexfiberをgeneric kernelの単純代入で代表させていない。**
この橋はprovenconstant-rank theoremからのvectorbundleであり、未証明flatnessを仮定しない。

λ無限遠は旧全base支持定理のW二sourceとして別に戻される。次節で全vertexを処理する。

## 4. W_ss/W_nilはEuler因子からprojective rank≤2

E=Σxi∂i、
`N=x0∂1+2x1∂2+3x2∂3`, `H=x1∂1+2x2∂2+3x3∂3`。
保存canonicalbasevectorを元九basisへ戻すと

`W_ss: P=E∧H`, `W_nil: P=E∧N`。

P=E∧Y と書き、E5=E+c∂cとすると、全vertex十八のままで

`Π=E5∧Y + (c*(Y+d)+V)∧∂c`。

projectiveEuler項を除いたtensorはdecomposableなのでrank≤2、actualDは恒等零。
Jacobi subsetやkernel四へ制限する前の全sourceに成立する。
無限遠をgenericλへ無理に代入した結論ではない。

## 5. B_nil全λ：actualquinticはhighest-weight productしかない

保存B_nil canonicalbaseは全λに対して、元全六括弧上で

`P=N∧Z, Z=−λE−H/2, d=(8λ+5)N/2`

である。したがって全canonicaltensorは

`Π=N∧(Z+(8λ+5)c∂c/2)+V∧∂c`。

signedEulerPfaffianはcに非依存である：Π²のc-dependent項は∂cを二度含んで零、
Eulerのc∂c項も同じ理由で消える。
canonical divergence-zero Poisson liftの一般Pfaffian–Casimir identity
`Π#dF=0`（既受理一般lemma）を使うと、

`N*Z(F)−Z*N(F)=0`。

λ≠0ではN/Zの0,1成分det=λx0²、λ=0では1,2成分det=−q0なので、
両linearfieldは全λでrational function field上独立。
従ってN(F)=Z(F)=0。
Fはdegree五homogeneousなので、`H(F)=−10λF`。

### Degree五のker Nを全て列挙する有限SL2証明

次を置く。

```
I=x0*x2−x1²,
J=x0²*x3−3*x0*x1*x2+2*x1³,
Delta=(J²+4I³)/x0²
     =x0²*x3²−6*x0*x1*x2*x3+4*x0*x2³+4*x1³*x3−3*x1²*x2².
```

Nはx0,I,J,Deltaを全て零にする。
degree五polynomialの次の六個はnonzero最高weight vectorsである。

`x0^5, x0³I, x0²J, x0I², IJ, x0Delta`。

Hのweightsは順に0,2,3,4,5,6。
SL2のCartan `3E−2H` のweightsは15,11,9,7,5,3。
対のloweringfield `3x1∂0+2x2∂1+x3∂2` とNは通常のsl2作用を与える。
各vectorが生成する非同型irreduciblemoduleのdimensionsは
16,12,10,8,6,4、その和56はdegree五四変数polynomial全空間のdimension C(8,3)に等しい。
従ってこれらが全分解を尽くし、ker Nの全六eigenlinesでもある。
finite highest-weight representationの標準sl2関係からの証明で、
未知Weitzenbock生成系や未実行kernel計算を仮定しない。

H(F)=−10λFだから、F非零ならこの六個のいずれかのconstant倍。
全てproper product又はpowerで、被約既約quinticは一つも無い。
他のλではF=0。
よって **B_nilの全λ・全vertex・全十九壁の当該補集合は対象空**。
特に旧kernel五のλ=−2/5、−3/5、−9/10も、順にx0I²、x0Delta、零として除かれる。
新M1や各wallのrankprofileは不要である。

## 6. なお残るexact whole sources

全sourceのactualnormalizer inverse、base零、W二源、Aの全U、B全lineを戻すと、
対象Eの未処理sourceは **A_ss/A_nil × λ2/5,3/5 の四全六kernel fiber**に限られる。
各sourceは保存済み元full350×18の完全kernelから全六coefficientを戻す。
この四源のactualD因子／named parent／空性を新しい有限whole評価で調べるのが次義務である。
次元上界だけで既知T38へ入るとはしない。新componentcountは本稿では零。
