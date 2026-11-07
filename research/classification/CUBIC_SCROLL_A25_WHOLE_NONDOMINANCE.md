# Smooth cubic scroll の残る A25 全源：quadric-net 境界による G 非優越

2026-11-07、G、producer plane_collision_review。
**純幾何証明。独立review前。** 新minor、GB、jet又はCASは実行していない。
固定scrollの射影Poisson像を六次元以下、全PGL5飽和像を二十四次元以下とする。
全点のE帰属、cubic cone、成分数又は全分類の完成は主張しない。

## 1. 全二型と、今回独立に読解した入口

S=F1=S(1,2)⊂P4、H=C0+2f、A=2C0+5fとする。
πは一般rank四、実Dは被約既約quintic、Sは非正規支持の一つ、π|Sは一般rank二。
R=div(π|S)の可約・非被約を許す。前段の全primitive root還元とcanonical-root排除の後に残るのは

* node型：G=γ ell1 ell2、二つのprimitive rows elli:N→A、div γ=C0。
* cusp型：G=γ ell0²、primitive ell0:N→A、div γ=C0。

である。これはnormal quadraticのrank二・一を指し、全Dの局所特異点が全域node/cuspとの仮定ではない。
Rは任意の重複成分を保持する。

今回、次の凍結producer notesを独立に全読解した。
`CUBIC_SCROLL_REMAINING_ROOT_SCALAR_CONNECTION.md`、
`CUBIC_SCROLL_QUADRIC_NET_GEOMETRY.md`、
`CUBIC_SCROLL_QUADRIC_NET_POISSON_REDUCTION.md`、
`CUBIC_SCROLL_EXCEPTIONAL_RESTRICTION_UNIQUENESS.md`。
以下に必要な全矢印を再記し、その先の新境界証明を与える。

rank二ではflat normal connection βに対し二root rowsのcharactersの差ηは、det rowsのprime multiplicity mと
Rのmultiplicity nからpole order≤n−mを持つ。
二rowsがprimitiveなのでcollision上のrank一matrixに対するprojectionのnilpotent off-diagonal entryはunitであり、
scalar characterによる相殺でこの上界を失わない。
R=C0+div det MからηはC0上の単極しか持てない。
closed性、residue class[C0]≠0、H0(Ω1_F1)=0によりη=0。

rank一ではE=P(N)上のclosed inverse residue αは唯一horizontal Γに高々二位、vertical C0に高々一位のpoleを持つ。
P1 fiberのresidue theoremでΓのresidue零、Gysin class[C0]でvertical residueも零。
Eはrationalなので、zero-residue closed rational one-formはrational primitive dfを持つ。
これはrational function fieldの一変数partial fractionsを順に積分し、必要なsplitting fieldからtraceする手証明である。
fのpoleはΓに高々一位だから

`f=ell1/ell0,  ell1∈H0(N*⊗A)`。

ここで**第二row ell1はprimitiveと仮定しない**。そのdivisorial zeroを全て許す。
dfのfiber成分が非零なので二rowsは独立。
inverse Pfaffianのvertical residue σ/(γ ell0²)とdfのfiber numeratorを比較すると
`det(ell1,ell0)=const*σ/γ`。
rank二でも同じdeterminant式が成り立つ。
従って両型で、二つの独立なglobal rowsとそのratio fを用いて

`σ=const*s_C0 det M,  ker α=ker df`

と書ける。nodeでは実際α=λdlog f、cuspではα=dfである。
nodeの任意GL2 row mixingは行わない。

## 2. Quadric-netの全P2-bundleと全full fibers

q=(x0x3−x1x2,x0x4−x1x3,x2x4−x3²)とし、net baseをB=P2_(a:b:c)とする。
Y=Bl_S P4は全incidence

`a x3−b x2+c x0=0,  a x4−b x3+c x1=0`

であり、全baseで係数rank二だからY=P(F)、
`0→F→O_B^5→O_B(1)^2→0`。
HはP4 hyperplane、Bもbase hyperplane classに用い、E=2H−B、D'=H+2Bである。
D'はF→O_B(2)というrelative hyperplane sectionで、全fiber P2となるzeroを許す。

c≠0ではplaneの座標X=x2,Y1=x3,Z=x4に対しEはQ=XZ−Y1²=0。
L=(c=0)では(x2,x3,x4)=t(a²,ab,b²)、Eは二つの異なるlinesに分裂する。
D'はV=C0により、そのfiberのC0 lineを含む。
s_D∈H0(F*(2))を五つのquadratics P_iで表示するとP0,P1はc倍。
二つのincidence syzygiesを一回引くことでP0=P1=0と出来、残る全表示は

`D': P(a,b,c) X+Q1(a,b,c) Y1+R1(a,b,c) Z=0`、P,Q1,R1は全quadratics。

この表示は一意であり、任意のextension coefficientを零と置いたのではない。
L上の全binary quarticは

`κ(a,b)=a²P(a,b,0)+ab Q1(a,b,0)+b²R1(a,b,0)`。

κ≡0ならp^{-1}LがD'の成分となり、Dの既約性に反する。
従ってκ≠0であり、長さ四の全zero schemeではfiberが丸ごとP2となる。
211/22/31/4の全partitions、L外の追加full fibersを許す。
normalizationのflat性はどこにも仮定しない。

## 3. 全net Poisson性とexceptional restrictionの一意性

一つのEuler gaugeで全raw Π#dq_iをI_Sに出来る。
その全域liftの理由はH1(T_S)=0とH0(T_P4)→H0(T_P4|S)→H0(N)の全射性である。
H1(T_S)=0はrelative tangent sequenceのO(2C0+f)のpushforward
O(1)⊕O⊕O(−1)から得る。
intrinsic conormal可換性をこのgaugeに戻すと{q_i,q_j}∈H0(I_S^(2)(4))。
全blowupで4H−2E=2B、p_*O_Y=O_Bなのでこの空間は正確にH0(O_B(2))の六次元。
従って全q-ratiosのbracketはbase quadratic bivectorを定め、p:Y→BはPoisson morphism。
base Poisson τが零なら二つのbase Hamiltoniansはgeneric P2 fiber上で独立なcommuting regular vector fieldsとなる。
両者はsmooth conicとlineを保つが、secant pairのLie algebraは一次元、tangent pairのBorelでも
commuting subspaceは一次元以下。rank四に矛盾するのでτ≠0。
Eのsingular conic fiberのnodeではdp|Eがrank一なのでτ|L=0。

さらにH0(Y,Λ²TY(−E))=0。
relative tangent filtrationのgraded piecesはO(H−B)、Trel(−2H+B)⊗p*TB、p*Λ²TB(−2H+B)。
後二つは各P2 fiberでT_P2(−2),O_P2(−2)のH0零、最初はH0(F*(−1))=0である。
最後はdual Steiner sequence O(−2)^2→O(−1)^5→F*(−1)から従う。
従ってπ_Eはπ_Y、元のπを一意に定める。
π_E自身はker dfとσで一意である。よって二rowsはprojective πを高々一つ定める。
scalar connectionや未読higher extensionを追加parameterとして数えない。

## 4. 全rowsの六係数表示とdual determinant

β:S→Z=P2のblowdown pointをp=(0:0:1)とする。
N*(2H)=β*Q_Z、A−2H=fにより

`H0(N*⊗A)=H0(Z,T_Z⊗I_p)`、dimension六。

全rowはprojective vector fieldのmatrixをscalarで割り、最後のdiagonal entryを零として一意に

```
T_i=[[A_i,0],[r_i,0]], A_i∈Mat2(C), r_i∈C²
```

と書ける。EはZのflag incidenceをSへ引き戻したもの。
generic net conicのparameter z=t/sを使うと

```
(U_i,V_i)=(a,b) A_i+c r_i,
l_i(z)=U_i z−V_i,
f=l_1/l_2 (node), または f=l_1/l_0 (cusp).
```

primitiveな元のrowに対してU_i,V_iにcommon linear factorは無い。
もし共通因子を持てば、固定conic parameterの零がruling fiber上のdivisorial zeroを与え、primitivityに反する。
rank一の第二rowにはこの主張を課さない。

Ω=U_1V_2−V_1U_2（cuspではindices1,0）と置く。
base Poissonの全式は非零constantを除いて

`τ=c Ω`。                                               (1)

次のdeterminant identityがその証明である。
Z上のrestriction cubicはF(y)=det(y,T_1y,T_2y)、dualはF*(b)=det(b,T_1^t b,T_2^t b)。
incidence b·y=0でf=(bT_1y)/(bT_2y)をannihilatorとするrank二bivectorを取り、ZへのpushをFとする。
そのdual baseへのpushはF*である。
例えばy=(x,y,1)、b=(a,−(1+ax)/y,1)でEの座標を(x,y,a)にすると

```
{x,y}=F, {y,a}=F*f_x/f_a, {a,x}=F*f_y/f_a,
F*(a*f_y−b*f_x)+y*f_a*F*=0.
```

最後の式は上の二つの三行列式を展開した恒等式であり、{a,b}=F*が従う。
generic affine chartで示した有理恒等式なので全coefficient・全境界へ延びる。
Tiの最後のcolumnが零だからF*=c det((U_1,V_1),(U_2,V_2))=cΩ。
τ≠0よりΩ≠0である。

両型のrelative hyperplaneは全表示として、constant倍を除いて

```
node: D'_fiber=V_1V_2 X−(U_1V_2+U_2V_1)Y1+U_1U_2 Z,
cusp: D'_fiber=V_0² X−2U_0V_0 Y1+U_0² Z.
```

Eのconic上でrootsが一致するのでこの式となる。
primitive rootのpair coefficientsのgcdは一であり、全coefficientsがdegree二なので余分なbase rational multiplierはconstantだけ。
κは、k_i(a,b)=bU_i(a,b,0)−aV_i(a,b,0)として、nodeではk1k2、cuspではk0²となる。
κ≠0よりnodeの両k_i、cuspのk0は非零。従ってscalar/zero A_iを黙って捨てていない。

## 5. Generic P2 fiber上の全connectionと非零constant

τがsymplecticなbase開では、Poisson submersionはflat horizontal connectionとparallel vertical Poissonに分解する。
horizontal vector fieldsはgeneric full P2 fiber上でregularであり、conic Eを保つ。
Aut(P2,conic)→Aut(conic)=PGL2は同型なので、conic上のf-horizontal distributionがこのconnectionを一意に決める。
linear matrixL=(U_i,−V_i)のSym²によるfiber座標変換でfをconstant conic coordinateにするとconnectionはtrivialとなる。
generic fiberのvertical Poissonは、nodeではstandard conic+secant line、cuspではconic+tangent lineを零因子とする。
Jacobiとτ≠0よりそのcoefficientはbase上の非零constant kである。
k=0なら全rankは二となりrank四に反する。

元の(X,Y1,Z)座標ではvertical cubicは、非零constantのnormalizationをkに吸収して

`(k/Ω) (XZ−Y1²) D'_fiber`。

これはGL2→Sym²GL2でconicがdet(L)²、projective volumeがdet(Sym²L)=det(L)³で変わることによる。
この表示はgeneric baseの全actual sourceに必要であり、全sourceをこの式から構成したとの逆主張は不要である。

ここで先に**c∤Ω**を全sourceで証明する。
Ω=c^m Ωtilde、m≥1と仮定する。
generic p^{-1}Lはactual E又はD'のdivisor成分ではないので、ωはそこでregularでなければならない。
§6の座標ではtrivial conic equationはq=c*v。
flat fiber conicのlocal defining functionのpullbackは、unitを除いてΩ²*qだからcのorderは2m+1。
κ≠0によりfのdenominatorはgeneric p^{-1}Lで非零、flat fiberのcollapsed conic pointは
standard D' lineの外にある（nodeでは両k_i非零、cuspではk0非零）。
従ってその近傍のstandard fiber two-formは(dQ0/Q0)∧regular-one-formの形であり、
pullbackはc上に高々logarithmic単極しか持たない。
一方base inverse formはτ=c^{m+1}Ωtildeなのでcにorder m+1≥2のpoleを持つ。
両者は相殺出来ず、ωのregular性に反する。
従ってc∤Ω、特にΩ_L≠0である。
これはgeneric coefficientへのspecializationの議論ではなく、全actual sourceのdivisor上のpole-order矛盾である。

## 6. 全L境界のclosed inverse residueが与えるbinary恒等式

§5でc∤Ωを全sourceについて証明した。
a=1 chartでLの座標をbとする。generic p^{-1}Lのdense chart x2=1では

`Y1=b+c x0, Z=b²+bc x0−c x1`,
`XZ−Y1²=c(−x1−b x0−c x0²)`。

従ってflat fiber conicのlocal equation Q0はgeneric p^{-1}Lでorder一で消える。
flat conic coordinate fのrestrictionは、nodeでk1/k2、cuspでk1/k0。
κ≠0によりstandard D' lineの外にあるgeneric pointなので、fiber two-formのresidueは
nodeでnonzero constant倍のdlog(k1/k2)、cuspでnonzero constant倍のd(k1/k0)である。
base inverse formのresidueは、τ=cΩのnormalizationでdb/Ω_L。
actual E+D'はp^{-1}Lをdivisor成分として持たないので、ωはそのgeneric pointでregular。
二つのresiduesを相殺すると、fixed nonzero constantsをkへ吸収して

```
node: Ω_L W(k1,k2)=k k1k2,
cusp: Ω_L W(k1,k0)=k k0²,       k∈C*.
```

Wはbinary quadraticsのJacobian（affine chartではk_i'k_j−k_i k_j'、固定数倍）。
両辺はhomogeneous degree四であり、generic Lでのrational differential恒等式から
**全L上のpolynomial恒等式**が従う。
Ω_Lの有限zero、κの全重複zero、根∞を除外したままにしていない。
これはparameter familyのgeneric係数から特殊係数へ推測したものでも、normalization-flat性からの帰結でもない。
各actual sourceのgeneric divisor上のclosed two-formを使い、恒等式を全curveへ戻した証明である。

## 7. Nodeの全binary quartic分類：22又は31だけ

k1,k2は非零binary quadraticsである。
比例ならW=0で右辺非零に反する。
また共通root pでは、(a,b)Aiが同じ(a,b)に比例するのでΩ_L(p)=0である。
**Wのdivisibilityだけでなく、このdeterminantの追加零を必ず使う。**

共通rootが無い場合、k_iのroot multiplicity mに対しord W=m−1だからord Ω_L=1。
Ω_Lはdegree二なので、κのdistinct rootsは高々二。
従って各k_iはdouble rootを持ち、その二点は異なる（κ partition22）。

共通rootが両方simpleならord W≥2、ord Ω_L≥1で左のorder≥3、右のorder二となり矛盾。
両方double又は二root共有なら比例してしまう。
残りは一方double、他simpleの共有でκ partition31だけ。
根∞もhomogeneous valuationで同じである。

traceを含む全matrixを戻す。
22型では座標を選び

`A1=λI+μ[[0,0],[1,0]], A2=ηI+ν[[0,1],[0,0]]`, μν≠0。

このときk1=μ b²、k2=−ν、
`Ω_L=λν+μν b+μη b²`。
binary恒等式から定数項と二次項が零、従ってλ=η=0。
両Aiは純nonzero nilpotent rank一である。

31型では共有rootをb=0に取り

`A1=[[λ,0],[μ,λ]], A2=[[a,0],[b0,d]]`, μ(a−d)≠0。

k1=μ b²、k2=(a−d)b+b0 b²、
`Ω_L=λ(d−a)b+(μd−λ b0)b²`。
恒等式の係数比較はλ≠0、d=0を強制する。
従ってA1はinvertible Jordan、A2はrank一semisimple。

両型でΩ_Lは二つの相異なる単根を持つ。
従って**この時点で**c∤Ω、Ωはglobalにもreducedと証明された。
double-line Ωやcを因子に持つ場合をresidueの段階でreducedと仮定していない。

## 8. Ωのinverse-form residueと31型の全排除

c≠0ではYはtrivial P2-bundleで、Ω-divisorはactual E+D'の成分ではない。
従ってinverse symplectic form ωはgeneric Ω×P2（E,D'の外）でregularである。
§5のflat座標ではωはbase inverse formとconstant standard fiber formの和。
standard conicのlocal defining functionは元座標でΩ²(XZ−Y1²)のunit倍となる。
そのためgeneric reduced Ω上でfiber formのresidueは、nonzero constantを除いて

`dlog f0`、f0はrank一matrixL|Ωのrow ratio。

より明示的にはstandard conic近傍でfiber formは(dQ0/Q0)∧α0+regular、
α0はnonzero constant倍のdlog f。
Q0のpullbackはorder二で消えるので、resΩは2α0(f0)となる。
base inverse formのresidueη_Bと相殺して

`η_B=nonzero_constant*dlog f0`。                       (2)

τ=cΩだからη_BはΩの各componentのgeneric点で非零で、polesはL∩Ωだけである。

Ωが二linesの和なら、linear map P2_B→P(Mat2)=P3のrankを分ける。
rank二以下では各Ω componentのimageはrank一の一つのpointで、f0はconstant。
rank三ならimage planeのSegre quadricとの交わりは二つのopposite ruling lines。
一方のrulingではrow ratioがconstantとなる。
どちらも(2)の左辺非零に反する。従ってΩはsmooth irreducible conic。

この場合Lのlinear matrix mapはplane embeddingであり、ΩはSegreの(1,1)曲線。
row ratio f0:Ω→P1はdegree一の同型である。
そのzeroとpoleは、それぞれ二row (U_i,V_i)が同時に零となるbase pointで、互いに異なる。
(2)からこれら二点はL∩Ωに含まれる。

しかし31型のinvertible A1では、(a,b)A1=0となるL上のpointは無い。
二linear forms U1,V1はindependent（primitive root）なので唯一のzero pointを持ち、それは必ずc≠0。
これは(2)に反する。
よって**31型の全sourceは空**、nodeに残る全Aiは純nilpotent rank一の二つだけである。
L∩ΩのtangencyやΩのnonreduced性を一般位置で捨てたのではなく、前節と(2)で全排除した。

## 9. Cusp：第二rowをprimitiveにせず全block五次元

ここではΩに重複lineやc因子があるかを事前仮定しない。§6のbinary恒等式だけを使う。
k0は非零。k1=0又はk1比例k0ならW=0となり不可能。
k0のsimple root pでk1も零なら、共通fixed-rootからord Ω_L≥1、ord W≥2となり右辺order二に反する。
k1非零ならord W=0なのでord Ω_L=2。
従ってk0が二つのsimple rootsを持てばΩ_Lのdegreeが四以上必要となり不可能。

k0はdouble rootを持つ。その点でk1が非零ならord W=1からord Ω_L=3が必要となり不可能。
k1もdoubleなら比例するので不可能。
従ってk0のdouble rootはk1のsimple rootであり、ord Ω_L=2である。

その点をb=0として全tracesを戻すと

`A0=[[λ,0],[μ,λ]], A1=[[a,0],[b0,d]]`, μ(a−d)≠0。

Ω_Lの一次係数は±λ(a−d)、二次係数は±(λ b0−μd)。
order二からλ=0、恒等式の非零constantからd(a−d)≠0。
従ってA0は純nonzero nilpotent rank一、A1はその固定pointを保つ全三次元linear spaceの中にある。
追加のopen conditionsを落とす必要は無く、上界にはその全三次元spaceを用いればよい。
global Ωのrepeated component、第二rowのdivisorial zero、a=0等もこの上界に含まれる。

## 10. 真のpresentation fibersと全G上界

nodeでは各nonzero nilpotent Mat2 coneのdimensionは二、r_iは各二係数。
全固定S pair source≤2+2+2+2=八次元。
二rootのindependent C* scalesはker dfを保ち、det Mとσをscalar倍するだけで、projective π_Eを保つ。
§3のglobal injectivityからprojective πも保つ。
二rows非零なのでこの二次元actionはfree、従って固定Sの射影像≤八−二=六。
nodeで任意GL2 row mixingを商にしていない。

cuspではA0 nilpotent二、A1のfixed-point-preserving space三、r0,r1各二で全pair≤九次元。
`ell0→a ell0, ell1→b ell1+c ell0`, ab≠0はker dfを保ち、σをscalar倍する。
二rows独立なのでこの三次元actionはfree。
第二rowがprimitiveでないsourceもこのactionで保持され、projective π像の真のfiberに入る。
従って固定Sの射影像≤九−三=六。
normalizationのflat性、generic四full fibers、sourceのsmoothnessを仮定しない。

最後にAut(S⊂P4)のdimensionは六。
PGL5×全固定source上の対角Aut作用は第一因子でfreeなので、source stabilizerが増えてもfiber六を失わない。
両型の全PGL飽和像は

`24+6−6=24`

以下となる。
これで残A25の全二型はGの未知dimension25–27成分の一般元から除かれる。
既済canonical-root branchesと合わせてsmooth cubic scroll rank二の全G非優越が従う。
全点の成分帰属、cubic cone、他の支持次数は別義務である。

## 11. 保存とreview scope

全κ partitionsを§6–9のhomogeneous identityとvaluationで処理し、roots at infinityも含めた。
全full P2 fibersはglobal smooth Yの計算内に残し、D'又はそのnormalizationをflatとした箇所は無い。
node Ω reducedは§7で証明した後だけ§8でresidueを用い、cusp Ω reducedは使わない。
原sourceの全extensionはexceptional restrictionのglobal H0 vanishingで処理した。

特に独立reviewでは、§5のc|Ω補集合、§4のdual determinant、§5のgeneric flat normalization、
§8のinverse-form residue、および§10の実presentation fiberを監査する。
数学的proofの成立とjob gate/hash一致は別であり、producer完成時点ではreview前である。
