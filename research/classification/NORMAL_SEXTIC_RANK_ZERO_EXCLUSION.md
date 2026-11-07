# normal degree六 rank零曲面の全排除：nef数値式とrelative conic

2026-11-07。主義務G、producer は plane_collision_review。
本稿は独立確認済みの `RANK_ZERO_NORMAL_SEXTIC_REDUCTION.md` を出発点とする新しい幾何証明。
rootによる別読解・受理前はproducer証明として扱う。
**命題。** normal・integral・非退化degree六曲面S⊂P4上で零となり、
実Dが既約被約quinticとなるrank四Poisson tensorは存在しない。
node、clean、normalization/conductor smoothness又はvertical divisor不存在を仮定しない。

## 1. 既受理の入口と下降

既受理の部分還元から、q(S)>0、sectional genus g≤3、全特異点rational、
minimal resolutionρ:X→Sは正種数b=q(S)の曲線B上のruled surfaceのblowupである。
smooth Sは既受理の小次数smooth irregular曲面補題でdegree五楕円scrollとなるため、
degree六の場合は既に除かれる。従ってSをsingularと仮定する。

`NORMAL_SEXTIC_RULED_NUMERICS.md` §1のproper graph論証によりrulingはg_S:S→Bへ下降する。
S→Bはflat、各fiberはCMの純一次元で、H=O_S(1)に対するHilbert polynomialはa*t+1。
a=1なら各fiberはscheme P1、S→Bがsmoothとなり矛盾するのでa≥2。
Riemann–Hurwitzとhyperplane完全列からg=2又は3となる。

## 2. elementary transformを仮定しない全nef数値式

L=ρ*Hとし、ε:X→X0をminimal ruled modelへのblowdownとする。
total exceptional basisで `L=ε*L0−Σm_iE_i` と書くと、nefnessとeffectivityから
`0≤m_i≤a`。infinitely near中心を含める。
ruled intersectionとprojection formulaから

`6=a(8−2g)+2a²(b−1)+Σm_i(a−m_i)`。

全summandが非負なので、b>1又はg=2は不可能。
残りはb=1,g=3であり、

`a=3, 全m_i∈{0,3}`、又は
`a=2, m_i=1が二個、他は0又は2`

だけである。詳しい式と全整数支持は前掲numerics noteに分離した。
Sへsmooth double-point formulaを適用していない。

## 3. a=3では全fiberがsmoothとなる

各vertical irreducible curveのnumerical classは、ruling fiber Fとtotal E_iの整数結合である。
L·F=3、L·E_i=m_i∈3Zなので、そのL-degreeは三の倍数。
Lはnefでtotal fiberのdegreeは三だから、各fiberにはL-positive componentが厳密に一つあり、
そのmultiplicityは一、L-degreeは三である。
他のcomponentはL-degree零なのでρによって点へcontractされる。

従ってSの各fiberのfundamental cycleは一つのirreducible curve、multiplicity一となる。
fiberはCMなのでgeneric reduced性からreducedであり、integral curveとなる。
flat Hilbert polynomialからarithmetic genusは零。integral projective curveでは
`p_a=g(normalization)+delta` なので、全fiberはsmooth P1となる。
flat morphismの全fiberがsmoothならg_Sはsmoothであり、Sもsmooth。
これは§1に反する。m_i=0/aをelementary transformで消す手順は不要である。

## 4. a=2の全fiberは、nonreducedの場合もplane conic

以下g=g_Sと略記する。ρのrational singularitiesと、p=gρがP1-bundleのblowupであることから

`ρ_*O_X=O_S, R1ρ_*O_X=0, R1p_*O_X=0`。

Lerayより `R1g_*O_S=0`。curve baseへのcohomology/base changeから、
**全fiber Fについて**h1(O_F)=0となる。geometrically reduced fiberだけに限定していない。

H_Fは五つのambient sectionでglobally generated。
`0→K→O_F^5→H_F→0` でH1(O_F)=0、dimF=1なのでH2(K)=0。
従ってH1(H_F)=0。Hilbert polynomial2*t+1からh0(H_F)=3である。
五つのrestriction sectionのspanが二以下なら、既にあるclosed embedding F⊂P4は
FをP1内へ埋め込むことになる。しかしP1内の純一次元closed subschemeはP1自身だけでdegree一。
これはdegree二に反する。従って五sectionは全三次元H0(H_F)を張る。

完全H_FでFはP2へclosed embedされる。CMの純codimension一schemeなのでCartier divisor、
degree二よりplane conicである。二直線又はdouble lineを先にreducedと仮定したわけではない。
この証明はribbonの外部分類を必要としない。

## 5. rank三bundleと負次数sectionによる矛盾

E=g_*Hと置く。前節の全fiber消滅とbase changeによりEはrank三vector bundle、
R1g_*H=0。さらに全fiberで五sectionが全H0を張るので、
`C5⊗O_B→E` はsurjective。従ってEはglobally generatedで、h0(E)=h0(H)=5。

quotient conventionのP(E)を用い、ξ=O_P(E)(1)とする。
相対evaluationと上の全section surjectionによりSはP(E)へclosed embedされる。
具体的にはP(E)⊂B×P4であり、元のSのgraphはこのplane bundleを通る。
P(E)はsmooth threefold、Sはintegral codimension一なのでCartier divisorである。
各fiberはconicだから

`S∼2ξ+π*D`, `H=ξ|S`

となる。e=degE、d=degDと書く。
ξ³=e、ξ²π*D=dより

`2e+d=H²=6`。

Bはellipticなので `K_P(E)=−3ξ+π*(detE)`。
Cartier adjunctionから

`K_S·H=−6+2(e+d)=d`。

g=3より左辺は−2。従ってd=−2、e=4。
elliptic Riemann–Rochはχ(E)=e=4を与え、h0(E)=5なのでh1(E)=1。
Serre dualityから非零map `E→O_B` がある。
そのimageはO_B(−Z)で、Eがglobally generatedなのでimageもglobally generated。
effective Zに対してこれはZ=0を強制し、E→O_Bはsurjectiveとなる。

このquotientが定めるsectionσ:B→P(E)上ではξはtrivial。
Sのdefining sectionをσへ制限するとH0(B,O_B(D))に入るが、degD=−2なので零。
従ってσ(B)⊂S。しかし `H|σ=ξ|σ=O_B` であり、Hのample性に反する。
これでa=2も排除される。

## 6. 到達点と未使用の文献

全normal degree六のrank零・既約被約quintic枝は、上記入口と二つのfiber次数の場合分けで空となる。
Halanayの1989年degree六分類本文、未取得publisher PDF、Reider、
elliptic rank二bundle分類、又はsingular surfaceへのsmooth double-point formulaは使っていない。
外部アクセス制限の再試行又は迂回もしていない。

これは非normal degree六、rank零plane、quadric cone、rank二restriction、又は全G/Eを尽くす命題ではない。
数値恒等式と有限支持の機械証明は別reportに限定し、fiber・descent・relative embeddingの
幾何学的な矢印を算術PASSで代用しない。rootの独立読解による受理後に限ってglobal主張へ接続する。
