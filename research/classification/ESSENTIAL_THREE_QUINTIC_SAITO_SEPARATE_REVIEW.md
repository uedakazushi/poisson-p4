# 本質三変数既約五次の全source・Saito上界：独立review

2026-11-07。主義務G。producerと別作者による数学的査読。**PASS**。

一般rank四Poisson構造の**実退化因子が既約被約quinticで本質三変数**なら、全canonical sourceはbalanced free plane quinticのcotangent表示と全六磁気項で尽くされ、PGL商前のprojective像は次元≤23である。
さらに、rank二planeに沿う位数四のleading symbolが`G₄=L*g₄`、g₄がsquarefreeなら、**実Dが既約quinticである場合は存在しない**。
後者はI,J非零性を要しない。特殊jのrank二pencil、重根、m2/m3、他のD支持及びE帰属は残す。

## 1. canonical liftを仮定で済ませない

五変数の二次bivector lift Pの自由度はE∧X（Xは一次場）である。
第一添字収縮をdivergenceとすると、div(E∧X)=5X−(tr X)E。
div Pは一次場で、div²=0によりtrace零なので、X=−div P/5でdivergence零代表Πが存在する。
核はscalar EulerだけでE∧E=0だから、この代表は一意である。

射影Poisson条件はE∧[Π,Π]=0。
BV恒等式とdiv Π=0からdiv[Π,Π]=0も成り立つ。
β=ι_[Π,Π]Ωは係数三次・形式次数二、Euler重み五で、ιEβ=dβ=0。
Cartan公式により5β=0、従ってΠは**strict Poisson**である。

α=ι_(Π²/2)Ωは係数四次・一形式、Euler重み五。
strict Jacobi及びdiv Π=0からdα=0なので、実五次式D=ιEαに対し

`5α=dD`, `Π#dD=0`。

最後は奇数次反対称行列のPfaffian核恒等式である。
これで任意のEuler代表をlogarithmic又はCasimirと置かずに、必要なcanonical恒等式を得る。

## 2. primitive gradientが全base括弧を消す

座標を(a₀,a₁,a₂,u,v)、D=F(a)とする。
Πのbase3×3 skew block AはA∇F=0を満たす。
符号付き三vector bで書けばb×∇F=0なので、fraction field上でb=c∇F。

Fが被約なのでgcd(F₀,F₁,F₂)=1。
共通primeはEuler公式からFを割るが、そのreduced factorのgeneric smooth pointでは全微分を割れない。
このprimitive性はS[u,v]への拡大でも保たれる。
各primeのvaluationを比較するとcの全valuationは非負、従ってUFD上cはpolynomialである。
bは次数二、∇Fは次数四なのでb=0。
ここではvertex依存を最初から除いていない。

よって全Πは

`Π=U∧∂u+V∧∂v+B∂u∧∂v`, `U(F)=V(F)=0`

となる。U,V,Bはこの段階では全五変数のquadraticである。
実Pfaffian/Euler恒等式は`D=−det(E_base,U,V)`。
u=v=0への評価から、baseだけのU₀,V₀について`det(E_base,U₀,V₀)=−F`を得る。

## 3. determinant criterionを代数的に自足証明する

ここでは取得出来ていないSaito1980本文に依存せず、必要なcriterionを直接示す。
M=(E_base,U₀,V₀)のdeterminantは−Fで、三columnはlogarithmic derivationである。
任意のδ∈Der(−log F)に対し、Cramer係数の分母はFだけである。
各reduced factorのgeneric smooth pointで、M mod factorはrank二。
rank一以下ならdet Mはそのfactorで二回以上割れ、F被約性に反する。
三columnは接平面に入りrank二なので、その接平面を張る。δも接平面に入るため、
Cramerの全numeratorは当該factorで割れる。
全factorについて割れば係数はSのpolynomialとなる。

従って`Der(−log F)=S E_base ⊕ S U₀ ⊕ S V₀`。
EulerはFへ5Fを作用させるため、そのkernelは厳密に

`AR(F)=S U₀ ⊕ S V₀`。

両基底は次数二なので、次数零・一次syzygyは存在しない。
U,Vをvertex monomialで展開すると、正のvertex次数の各coefficientはbase次数一以下のgradient syzygyとなり、全て零。
したがってU=U₀,V=V₀であり、全fiber依存を後から証明して消した。

## 4. Jacobi・divergence・全磁気項

[U,V]=aU+bVは次数比較からa,bがbase一次式である。
`det(E,U,V)=−F`のLie微分をU,Vで取ると、a=div V、b=−div U。
strict Jacobiの(aᵢ,u,v)成分は

`[U,V]+B_u U+B_v V=0`

なので、独立性より

`B=−(div V)u+(div U)v+R(a)`。

Rは**任意のbase quadratic六係数**である。
残るJacobiはbase bracket零及びU,Vのvertex非依存から自動零。
このBは同時にdiv Π=0を満たし、必要条件だけの選択sourceではなく完全な表示になる。

固定projective Fの全U,VはAR₂の全frame GL₂、RはA⁶である。
同時tensor scaleを除くとfiber dimensionは4+6−1=9以下。
frame determinantによるDのscalarはprojective Fに含まれ、別parameterではない。

## 5. Jacobian schemeはperfect・saturatedでτ=12

AR基底から、J_F=(F₀,F₁,F₂)に対するgraded resolutionは

`0 → S(−6)² → S(−4)³ → J_F → 0`。

S/J_FのHilbert numeratorは1−3t⁴+2t⁶で、t=1の零位数二、二階微分24。
したがってKrull次元一、Hilbert多項式は12。
同じminimal free resolutionからdepth一でCohen–Macaulay。
irrelevant idealに支えられるtorsionがなく、J_Fはsaturatedである。
FはEulerによりJ_Fへ入るので、射影有限scheme Σは各特異点のTjurina schemeそのもの。
よって全τ(F)=12。

五twist後は

`0 → O(−1)² → O(1)³ → I_Σ(5) → 0`。

P²上H¹O(1)=H²O(−1)=0よりH¹I_Σ(5)=0、h⁰I_Σ(5)=9。
曲線F自身もこのspaceに入り、projectivization後のequianalytic tangent上界は八である。

## 6. µ=τと有限analytic型：正確な外部入力

F既約quinticの算術種数は六。
Milnor–δ公式から、gを正規化種数、r_pを枝数として

`µ(F)=12−2g−Σ(r_p−1)≤12`。

各点でτ_p≤µ_p、全τ=12なので、g=0、全r_p=1、各点でµ_p=τ_pとなる。
孤立germの有限Milnor代数M=O/(f_x,f_y)とM/(f)の長さが等しいため、f∈(f_x,f_y)。

ここで使うSaito1971の標準定理は次の**正確なscope**である。
孤立特異点を持つconvergent complex hypersurface germ fについて、f∈J(f)なら、解析座標変更後fはweighted homogeneousになる。
独立にMaynadier仏訳のindexed Introductionから、この仮定と結論を読解した。
これはfree divisor determinant criterionとは別の定理であり、全文証明を今回再証明したとはしない。
元の独語本文は未取得、仏訳全文openは11,602,116 bytesのreader上限で失敗、1980原文openは429であった。
これらを`CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_PRIMARY_SCOPE_20261107.json`に区別して保存した。

weighted homogeneous・irreducible・unibranch plane germは、unit及び解析座標を除きxᵃ+yᵇ、gcd(a,b)=1、a,b≥2である。
実際、coprime weights上の同次数monomialは一つのlattice segmentに乗り、
coordinate monomialを除けばbinomial factorsの積になる。irreducibilityは一factorだけを許す。
そのµは(a−1)(b−1)。µ≤12から、順序を除く全有限supersetは

`(2,3),(2,5),(2,7),(2,9),(2,11),(2,13),(3,4),(3,5),(3,7),(4,5)`。

このlistの全点がquintic上で実現するとは主張しない。有限性に必要なsupersetとして使う。
全µ=12を持つ局所型multisetの候補は十九通りであり、点数は六以下。

## 7. 有限strataと「位置移動を含む」接空間

上の有限個のisolated analytic型には共通の十分大きいcontact-determinacy次数kを選べる。
さらにkを大きくし、各型でm^(k+1)⊂(f,f_x,f_y)となるよう取る。
有限k-jetのcoordinate-change/unit群は代数群で、その各orbitはlocally closed。
したがって、特異点の位置もincidence parameterとして動かし、所定のjet-orbitに入る条件を課せば、
全F族は有限個のconstructible equianalytic strataに入る。
点数と型multisetは有限で、追加特異点の不存在もconstructibleに扱える。
ここでは全rational cuspidal quinticの分類や、PGL軌道の有限性を仮定しない。

各stratumをsmooth locally closed strataへ細分し、その一般点でtangentを取る。
jet-orbitの曲線は局所的に群へlift出来るので、各特異点でのFのfirst derivativeは
(f,f_x,f_y)に入る。点の移動はf_x,f_yへのconstant vector field項を与える。
**点を固定した条件のmJではなく、位置移動を含むTjurina ideal Jを使う**。
高jetの残りも選んだm^(k+1)包含で同じidealへ入る。
従ってembedded projective tangentはH⁰I_Σ(5)/C Fへ入り、dimension≤8。
有限stratification全体にもこの上界が及ぶ。

Fを動かす際、AR₂はlinear syzygy coefficient matrixのrankが一定十六のkernel bundleで、rank二。
そのminor chartsでalgebraic frameを選べるので、基底選択に隠れた追加dimensionはない。
本質三変数WのchoiceはGr(3,5)の六。
各chartで固定した補空間に対して全Πを既に表示したので、別にvertex shear六を加えない。
実際u'=u+ℓ(a),v'=v+m(a)はU,Vを保ち、

`R' = R+(div V)ℓ−(div U)m−U(m)+V(ℓ)`

という別の任意quadraticへ変えるだけである。
よってPGL商前の全projective像は**6+8+9=23**以下となる。

## 8. ordinary四重点によるscalar m4全枝の排除

rank二plane沿いnormal order四で、leading G₄=L*g₄、g₄ squarefreeとする。
五次斉次性だけから**全**Dは

`D=L*g₄(z,w)+g₅(z,w)`

である。未知g₅を固定していない。
実Dが既約なら本質二変数以下は不可能なので、これは本質三変数の既約quintic。
p=[L:z:w]=[1:0:0]は四相異接線のordinary quadruple pointである。
δ_p=6、r_p=4、µ_p=9。
算術種数六をこの点が使い切り、他の特異点は無いので全τ≤全µ=9。
一方、本質三変数rank四Poissonの前節はτ=12を強制する。矛盾する。
同じ矛盾は§6の「全特異点unibranch」にr_p=4が反することからも得られる。

このargumentはg₄のI,J非零性を使わないため、**特殊jのrank一leadingも排除する**。
既受理parallel lemmaがG₄=L*g₄を強制するconcurrent・genericjの旧上界31枝も、D既約なら空となる。
新たなminor、P1/P2解、追加jet又はGBは不要である。
特殊jのrank二pencilは全Dが本質四変数となり得るので、この矢印には含めない。

## 9. 証明書と検証の範囲

独立proof auditを本文に記し、finite job gate後の標準ライブラリ整数計算は、
Hilbert numerator/τ12、h⁰I5=9、genus/branch budget、全十binomial型と十九multiset、
6+8+9=23、ordinary4pointのδ6/µ9だけを検証した。
これは外部解析定理を計算で証明したという意味ではない。

|artifact|SHA256|
|---|---|
|producer ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md|8a941b6a436a118d9833cd3cb912ca3e448a048c88ee4e386019e027da85a217|
|code/classification_essential_three_quintic_review.py|67ce38c5fb08261334ed91ce0bc0f75421ca778d62db60c29d72e032daccbc55|
|reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_REVIEW_20261107.json|350e5e07a0ad98dce37702cd1f25727f967a501bde3d4d6d6a50b2690f56a716|
|reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_REVIEW_JOB_20261107.json|bcfda47e7b05a8a31035df7d3c819aafce836cd56f82f71fd85582ec20bfa9f8|
|reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_PRIMARY_SCOPE_20261107.json|5754c56a8081c5eff6b33a14db30c103161f344c5d3b9f3f02282f81a8f13afa|

全依存は公開path/hashで固定した。global台帳・成分数・原稿・Gitは変更していない。
