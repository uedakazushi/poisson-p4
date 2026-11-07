# normal次数六rank零の全排除：統合担当による独立幾何確認

2026年11月7日。producer は plane_collision_review、reviewer は root。
`NORMAL_SEXTIC_RULED_NUMERICS.md` と `NORMAL_SEXTIC_RANK_ZERO_EXCLUSION.md`
の全証明を別に読解し、**全normal・integral・非退化degree六rank零枝の排除として受理する**。
出発点の部分還元は既受理の `RANK_ZERO_NORMAL_SEXTIC_REDUCTION.md` とその別作者reviewである。
新たな文献分類、clean、node又はvertical不存在の仮定は付けない。

数値恒等式を独立に導くと、minimal ruled modelへの全total exceptional basisについて

\[
 L_0^2=6+\sum m_i^2,\quad
 K_XL=-L_0^2/a+2a(b-1)+\sum m_i=2g-8
\]

なので `6=a(8-2g)+2a²(b-1)+Σm_i(a-m_i)`。
total exceptional divisorもtotal fiberとの差もeffectiveであり、nef Lから
`0≤m_i≤a` はinfinitely near中心にも成立する。
正種数baseへのrulingの下降、fiber degree一のsmoothness、hyperplane完全列と
Riemann–Hurwitzを合わせると、残るのは b=1,g=3 と二つの数値型だけである。

a=3では全vertical degreeが三の倍数。各S-fiberは唯一multiplicity一の曲線となり、
CMとgeneric reduced性からintegral、flat arithmetic genus零からsmooth P1である。
よってS自体がsmoothとなり、既受理の小次数smooth irregular曲面定理に矛盾する。

a=2ではrational resolutionとP1-bundleのblowupから `R1g_*O_S=0`。
flat curve familyのtop cohomology base changeで全fiberのH1(O)が零となる。
五生成sectionのkernel完全列からH1(H_F)=0、したがってh0(H_F)=3。
元のembeddingはこの三sectionを全て使う。rank二以下ならdegree二の純一次元schemeを
P1へclosed embedすることになり不可能である。従ってreducible/nonreduced fiberも
全てschemeとしてplane conicとなり、reduced性の追加仮定は不要である。

quotient conventionで `E=g_*H` はrank三globally generated bundle、h0(E)=5。
全fiber evaluationからgraphはP(E)のclosed subschemeであり、smooth threefold内の
integral divisorとして `S=2ξ+π*D` はCartierとなる。
`2degE+degD=6` とadjunction `K_SH=degD=-2` からdegE=4。
elliptic RRとSerre dualityから非零 E→O_B が存在する。globally generatedなimageは
O_B(-Z)なのでZ=0、したがってこれはquotientである。
対応するsection上でSの式は負次数Dのsectionとなり零、ゆえにsectionはSに含まれる。
しかしHの制限はtrivialとなりample性に反する。全てのfiberと全境界を保持した矛盾である。

有限算術reportは数値式と整数支持だけの検証として扱い、この幾何証明の代用にはしない。
Halanayの未読分類、Reider、rank二bundleの分類及びsingular surfaceへのsmooth
double-point formulaを使わないことも確認した。
非normal degree六、quadric cone、plane又はrank二restrictionは残る。
