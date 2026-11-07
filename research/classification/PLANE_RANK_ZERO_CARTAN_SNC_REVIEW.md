# rank零planeのCartan first symbol：全higher源のMPS40帰属を独立確認

2026年11月7日。主義務G。**producer命題を受理する。**
rank零plane、m_D=2、非零semisimple constant normal sl2 symbol、
共通global field Vとtrace field Tが同時対角化可能、という全枝が対象である。
このsymbolを持つ全higher Poisson tensorと、それを通る全global irreducible componentはMPS40の和に含まれる。
generic eigenvalue、nonresonant period、higher-weight消去、局所成分一意性は仮定しない。

独立checkerは `code/classification_plane_rank_zero_cartan_snc_review.py`、
reportは `reports/CLASSIFICATION_PLANE_RANK_ZERO_CARTAN_SNC_REVIEW_20261107.json`、
jobは `jobs/plane_rank_zero_cartan_snc_review.json`。
別作者標準Fraction helperだけを用い、producer/CAS/SymPyをimportしない。
gate後0.0053秒、最大RSS12,032KiBで全記号恒等式PASS。

## 1. 全weight零条件と実SNC端点

constant semisimple normal matrixをGL2で対角化し、その非零scaleをVへ吸収する。
全weight零sourceは

    Π0=zX∧∂z+wY∧∂w+(az²+bzw+cw²)∂z∧∂w,
    X=T+V, Y=T−V.

checkerは一般tracefree X,Yの全十六係数、a,b,cの全三係数を保ち、strict Jacobiの全非零vector係数が

    zw[X,Y]+(az²−cw²)(X−Y)

に一致し、他の七成分が零であることを再構成した。
projective Jacobiはこの式とEuler_baseのwedgeである。
X−Yが非零projective fieldなので、z²,w²係数からa=c=0。
zw係数はcommutatorがscalar Eulerであることを意味するが、trace零によりそのscalarも零。
従って全sourceで[X,Y]=0、K=bzwである。
traceless symbol自体が零の場合はX=Yでleading Pfaffianが零となり、m_D=2の対象ではない。

同時対角化後、X,Yの全六固有値とbを記号として保持した。
全350 projective Jacobiは恒等零であり、五つの符号付きPfaffianのEuler収縮は正確に

    D0=−κ*x0*x1*x2*z*w,
    κ=det[[1,1,1],[α0,α1,α2],[β0,β1,β2]].

normal二次symbolが非零というm_D=2の条件はκ≠0そのものである。
五つのhyperplaneは独立であり、D0は被約SNC。
Π0はtoric logarithmic tangent bundleのconstant bivectorで、κ≠0はその非退化性である。
従って単にSNC式を持つだけでなく、全域でlog symplecticである。
b又は一部biresidueの零を反転したり、特別periodを捨てたりしていない。

## 2. 全higher源を運ぶ実PGL作用

rank零性とconormal abelian性によりnegative normal-weight係数は無い。
Euler商後の全quadratic sourceはΠ0+Π1+Π2であり、weight一三十、二九の全係数を持つ。
checkerはこの全39 basisと全150 homogeneous coefficientのnormal weightを別に列挙した。

    g_t=diag(1,1,1,t⁻¹,t⁻¹),
    (g_t)_*Π=Π0+tΠ1+t²Π2

というpushforward指数を全係数で確認した。
右辺は全tで多項式、t=1が元のtensor、t=0が上のSNC tensorである。
t≠0で実PGL同値であり、Schouten括弧の自然性からPoisson、
従ってpolynomial恒等式としてt=0でもPoissonである。
formal jet、未証明のflatness、higher-weightの消去を使わない。

## 3. 境界所属の逆推論をcomponent安定性で補う

元Πを通る任意のglobal irreducible component Cを取る。
connected algebraic group C*とCの積は既約である。
その作用像閉包は既約で、identityからCを含むため、componentの極大性からCそのものに等しい。
従ってCは全t≠0のΠtを含み、閉性によりΠ0も含む。

重要なのは、**元Πを通る同じcomponent C自身が実SNC点を含む**ということである。
別の既知成分が境界点Π0を含むことだけから、元Πの所属を推論していない。
Π0での成分一意性も不要であり、Πを通る全componentにこの議論を適用できる。

## 4. 外部定理の全period scopeを一次資料で再確認した

Matviichuk–Pym–Schedler, arXiv:2010.08692v2の
[Theorem1.4](https://arxiv.org/html/2010.08692v2#S1.Thmtheorem4)、
[§5.3](https://arxiv.org/html/2010.08692v2#S5.SS3)、
[Example6.24](https://arxiv.org/html/2010.08692v2#S6.Thmtheorem24)を2026年11月7日に独立に直接再読した。
[PDF](https://arxiv.org/pdf/2010.08692v2)では印刷7頁、39–41頁、55–56頁である。

Theorem1.4はP4のPoisson moduliのうちSNC log-symplectic locusと交わるcomponentを四十個とする。
Example6.24ではprojective toric boundaryについて、全periodでresonant local-system obstructionが消えることを示し、
n=2について四十familiesのglobal distinctnessを直接確認してTheorem1.4へ接続している。
従ってTheorem1.2の一般nonresonant openを、未確認の特殊periodへ拡張していない。
§5.3のtoric log formの非退化性も、特定entryではなく行列全体の可逆性である。

P4のrigidityによりpair deformationでunderlying varietyの新modulusは生じない。
固定P4のPoisson parameter spaceからPGL5 quotientへのsmooth presentationはconnected groupの既約fiberを持ち、
上記componentと固定空間のcomponentが対応する。
scalar projectivizationもconnected Gm作用なのでcomponent対応を変えない。
既存T01–T40との内部ID照合は、受理済み `MPS40_CROSSWALK_20261031.md` と対応reportを依存とする。
今回その四十図式の抽出又は外部定理全体を再証明したとはしない。

よって§3のCは四十componentの一つであり、元ΠもMPS40へ入る。
nilpotent constant normal matrix、commutingだが非同時対角化のT,V、他co-Higgs型、m_D≥3は補集合として残す。
非SNCの既知境界へさらに退化する場合に、この論理をそのまま適用してはいけない。
