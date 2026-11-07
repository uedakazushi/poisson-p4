# degree六rank零支持の全排除：quintic正規化のsectional genus零

2026-11-07、主義務G。producerによる新しい幾何証明。
独立review前は本noteの結論をglobal完了主張へ混ぜない。

**命題。** D⊂P4をrank四Poisson tensor πの実退化因子である既約被約quinticとする。
非正規支持にdegree六のintegral surface Sが現れ、π|S=0なら矛盾する。
Sのnormality、smoothness、ordinary-double性、vertical divisor不存在又はcone性を仮定しない。
既受理の全normal-sextic排除を再計算せず、D自身の正規化から全非normal Sも扱う。

## 1. 全normalization-lengthをdegree六支持が使い切る

ν:X→Dを正規化、L=ν*O_D(1)とする。
`IRREDUCIBLE_NONNORMAL_QUINTIC_ENTRANCE.md` §1の全必要条件は

`Σ_i δ_i degS_i=6−g`, `L³=5`, `K_X+C_X∼0`

を与える。δ_iは正規化商の曲面一般点での正のlength、g≥0は一般平面五次切断の
正規化のgenusである。
degS=6なのでg=0、δ_S=1であり、二次元の非正規支持はSただ一つである。
導手Weil divisor C_Xの全成分はSへfiniteに写る。

このgをpolarized pair (X,L)のsectional genusと同一視する点を明記する。
Xはnormal threefoldなのでSingXのdimensionは高々一。
νはfiniteなので一般P2⊂P4はν(SingX)を避ける。
pullbackした二つの一般hyperplaneはXのsmooth locus上でBertiniによりsmooth complete-intersection
curve Γを作る。νはfiniteでDのdense open上で同型なので、Γは平面quintic切断Cへfinite birationalである。
一般Cはintegralであり、Γに追加のcurve成分は無い。
XのS2性により二hyperplaneはregular sequenceで、ΓはCMかつgeneric reduced、従ってreduced。
よってΓはCの正規化である。
二回のCartier hyperplane完全列、又はHilbert polynomialの二階差分により
sectional genus g(X,L)はp_a(Γ)=0に等しい。
これはX全体のsmoothness又はCのsmoothnessを仮定した同一視ではない。

## 2. all-normal-scopeのsectional-genus零定理

直接読解した一次資料は Andreas Höring,
*The sectional genus of quasi-polarised varieties*,
<https://math.univ-cotedazur.fr/~hoering/articles/a12-sect-genus.pdf>、Theorem1.2、印刷2頁である。
同定理はnormal complex varietyとnef・big Cartier divisorの全場合に対し、
sectional genus零ならΔ-genus零とする。
smooth、terminal、klt又はCMという追加仮定はstatementに無い。
Lはfinite pullbackなのでampleであり、当然このscopeに入る。
従って

`h⁰(X,L)=dimX+L³=8`。

完全linear system |L|は元の五sectionを含むのでbasepoint freeであり、
X→P7のmorphismは元のνをlinear projectionとして因数分解する。
Lがampleなのでこのmorphismはfinite、νがbirationalなのでそれもbirationalである。
imageはdegree五、dimension三、codimension四のminimal-degree varietyである。
minimal-degree varietyはnormalなので、finite birational性からXはそのimageと同型。
従って|L|は実際にembeddingを与える。

既読のminimal-degree分類（Brodmann–Park arXiv:0808.0090v3,
Notation and Remark2.1、印刷4頁）により、degree五threefold Xは

`S(1,1,3), S(1,2,2), S(0,1,4), S(0,2,3), S(0,0,5)`

のいずれかである。
最初の二つはsmooth rational scroll、残る三つはpositive-dimensional vertexを持つconeである。
quadric、projective space、Veronese及びそのconeはdegreeが違うので入らない。

## 3. singular scroll normalizationならSはconeになる

ν:X⊂P7→D⊂P4は五sectionによるlinear projectionであり、finiteなのでcentreはXを避ける。
Xがvertex Vを持つconeなら、そのvertexは射影で同dimensionのlinear space V'へ写る。
joinをlinear projectionすることにより、DもV'をvertexとするconeになる。
非正規支持の唯一の二次元成分Sは、このconeの非正規支持なのでV'をvertexとするconeである。
具体的にはDのdefining polynomialをvertex方向に依存しないものとして書けば、
そのnormalization-defect支持もbaseの非正規支持とのjoinとなる。

Sが非退化なら `RANK_ZERO_HIGH_DEGREE_CONE_EXCLUSION.md` の全source因子証明で既約quinticに矛盾。
Sがhyperplane内なら、degree≥3の全restriction-rank排除
（`IRREDUCIBLE_NONNORMAL_DOUBLE_SURFACE_DEGREE_RESTRICTIONS.md` のEuler/prime-CI-square補題）に矛盾する。
従ってsingular Xの全三型は除かれる。

## 4. smooth Xではδ一の二つの横断型を尽くす

残るXはP1上のsmooth rank三scrollで、PicXはLとfiber divisor Fで自由生成される。
standard projective-bundle canonical formulaから

`K_X=−3L+3F`, `C_X∼3L−3F`。

δ_S=1であるreduced plane-curve singularityは、generic algebraic closure上で
ordinary node又はordinary cuspである。
これはδ一ならmultiplicity二、node又はunibranchのsemigroup〈2,3〉だけ、という曲線局所計算である。
ordinaryを全域に仮定せず、degree六length条件から横断一般点の二型へ到達した。

generic cuspの場合、normalizationのconductor orderは二で、Sの上の唯一のprime divisor Γについて
`C_X=2Γ` となる。
XはsmoothなのでΓはCartier。従ってC_XのPic classは二で割れるはずだが、
`3L−3F` のL係数は奇数であり矛盾する。
一般点以外でcodimension二の導手schemeが変化しても、このWeil-divisor parityには影響しない。

## 5. generic nodeならPoisson bivectorが全導手で零となる

generic nodeの場合、C_Xはreduced Weil divisorである。
πはPoissonなので実退化因子Dのidealを保存し、Dのsmooth locusにregular bivector π_Dを誘導する。
そのrational pullbackをπ_Xとする。
大域的な「正規化上のPoisson bracket延長」定理を追加で仮定せず、全codimension一で直接確認する。

導手の一般点では、finite separable base change及びcompletionの後、Dは
`K[[u,v]]/(u*v)`（KはSのfunction fieldのseparable extension）となり、
各normalization branchはK[[u]]又はK[[v]]である。
π|S=0なので、D上の二つのregular functionのbracketはideal(u,v)へ入る。
branchのminimal primeはderivationで保存される。
実際δ(u*v)=0を(u)で割るとv*δu=0なのでδu∈(u)、他方も同様である。
各branchへのring mapはsurjectiveであり、その誘導bracketはbranch parameterで割れる。
従ってπ_Xは各reduced conductor branch上で零となる。
separable base changeでの零性はfaithfully descentする。

C_X以外のcodimension一ではνは同型でDはsmoothなのでπ_Xはregularである。
Xがsmoothであるため、全codimension一でのregular性とC_X上の一次零性は
Hartogsにより

`π_X∈H⁰(X,Λ²T_X(−C_X))`

を与える。
dimension三とK_X+C_X∼0から `Λ²T_X(−C_X)≅Ω¹_X`。
しかしsmooth rational scrollではH⁰(Ω¹_X)=0である。
一方、実DがreducedなのでDの一般点でπが零とはならない
（全係数がDのidealへ入ればPfaffianはその平方で割れる）。
従ってπ_Dはrank二、π_Xは非零。
これは矛盾する。

generic cuspとnodeの双方が除かれ、命題が従う。

## 6. 射程

本証明はdegree六rank零Sについてnormalと非normalの全場合を排除する。
normal-sextic noteのruled数値式、Halanayの未読分類、
surface Sの正規化・導手のsmoothness、又はconductorの全解析型分類を用いない。
degree五Sはlength残余とsectional genus一の可能性を持つため、この証明だけでは閉じない。
rank二restriction、plane、quadric surface及び全G/Eの残余は保持する。
