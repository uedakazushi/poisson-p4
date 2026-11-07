# 被約可約R・plane位数三：全root coversとprimitive degreeによる接続入口

2026-11-07、PRIMARY G、producer plane_collision_review。**純幾何の部分還元、独立review前**。
新source probe、minor、GB、jet又はCASは実行していない。
全G≤24又はE帰属はまだ結論しない。表の未処理枝を、既に終わったm二/m四又は非被約Rへ混ぜない。

## 1. 対象と全未処理表

πはP4のrank四Poisson、actualDは被約既約quintic、S=P2⊂SingDはrestriction一般rank二。
πSの零cubic Rは**被約可約**、S沿いnormal位数は正確に三とする。
N=O(1)⊗C²、normal variablesはu,v、

`G3∈H0(O_S(2))⊗Sym³ C²`, G3≠0。

Hamiltonian associated gradedによりγ=G3/fはflat conormal connectionのparallel tensorである。
ここでfはRのcubic equation。ordinary covector connection βはEuler β=I、fβの係数はdegree二polynomial。
R被約なので各entryのdivisorial poleは高々一位。
全matrixをlogarithmicとは仮定しない。

| R | 方程式 | π1(S\R) | SFで許すroot monodromy | 現在のm三義務 |
|---|---|---|---|---|
| transverse conic+line | z(xy−z²) | Z | 1,C2,C3 | SFの全covers、21全四degree flags、3 |
| tangent conic+line | z(xz−y²) | Z | 1,C2,C3 | 同左 |
| nonconcurrent three lines | xyz | Z² | 1,C2,C3 | 同左 |
| concurrent three lines | xy(x−y) | free group F2 | 1,C2,C3,S3 | 同左。S3 connected coverも保持 |

SFはgeneric partition111、21はdouble+simple、3はtriple rootを指す。
上表は全未処理branchのledgerであり、全monodromy型の実在を主張していない。
既約Rの全空性、m四の全R、被約R/m二SFのMPS被覆は既済で、本noteでは再開しない。
非被約R/m三はgeneric_exhaustion、被約R/m二doubleはquartic_sourcesが別担当である。

## 2. π1と全SF coversの有限性

tangent conicではz=1でU={(x,y):x−y²≠0}≅C*×C。
triangleはU≅(C*)²、concurrentはy=1でU≅(C\{0,1})×C。
transverse conicではz=1、U={xy≠1}。
x≠0の開は座標(x,xy−1)で(C*)²となる。
x=0のlineの小tubular近傍はlineへretractし、そのpunctured近傍はnormal meridian xを生成する。
このmeridianはxy−1=-1のままなので、van Kampenはxのgeneratorだけを殺し、残るπ1はZとなる。
lineの全有限yを含む近傍を取り、無記録の追加generatorを落としていない。

SFの場合root incidenceはU上finite étale degree三。
そのpermutation imageはS3のsubgroupで、全subgroup型は1,C2,C3,S3。
abelian π1の前三型ではS3 imageは無く、concurrentでは排除しない。
各固定Rでπ1は有限生成なので、S3へのhomomorphismsは有限個。
従ってsplitting Galois coversは有限個のmarked coverからなり、そのdegreesはN=1,2,3,6である。
connected degree三をC3だけと呼びS3を捨てるのは、concurrentの場合に誤りである。

## 3. Nontrivial SF coverでも使える全degree budget

上のsplitting coverをP2上でnormal compactificationし、さらに三root rational mapsを同時にresolveする。
得たsmooth surfaceをT、P2へのgenerically finite mapをp、H=p*O(1)とする。
各root map T→P1が与えるglobally generated line bundleをL1,L2,L3とする。
pullback G3は三root linear factorsの積に分かれ、normal DVRでcontentの積を取れば

`p*G3=h ell1 ell2 ell3`, `h∈H0(O_T(2H−L1−L2−L3))`。

hにdivisorial poleは無い。各primitive pairにはDVRでunit成分があるためである。
Hはnef、Liもnefなので

`sum_i Li·H ≤2H²=2N`。                                 (1)

これをroot degreesの全budgetとして用いてよい。
cyclic connected coverでは三値が等しく各≤2、S3 splitting degree六では各≤4、
1+2 cover degree二では固定rootの値＋moving二rootの同値の和が≤4。
ただしcompactificationのPicが零、全rootがpolynomial vector、又はnormalizationがflatとは仮定しない。
このbudgetだけでnontrivial cover上の全section familiesを分類したとは言わない。
exceptional curvesのH-degree零を捨てず、basepointsはresolution内に保持する。

## 4. Fully split SFの全四degree flags

rootsがC(S)上全てsplitすると、UFDと二gradingにより

`G3=h ell1 ell2 ell3`, `deg ell_i=d_i≥0`, `deg h=2−sum d_i≥0`

とglobalに書ける。各rowはprimitive homogeneous pairである。
constant normal GL2を使い、順序を除いて全degree flagsは

`(0,0,0), (0,0,1), (0,0,2), (0,1,1)`。                (2)

hと各nonzero pair determinantのdivisor支持はRに含まれる。
U外へpole/collisionが出ることはparallel nonzero stable tensorに反する。
各場合の全covariant入口は次である。

* 000：G3=h g0、g0はconstant SF cubic、hはR支持の全quadratic。
  三constant root linesを保つβはscalar、`β=dlog(f/h) I`。
* 001：二constant rootsをu,vとしG3=h uv(a u+bv)。
  a,bは相異なるRのline components、hもline component。従ってconic二型では存在しない。
  全βはdiagonal、entriesは`dlog(f/(ha)), dlog(f/(hb))`。
* 002：G3=uv(Au+Bv)、A,Bはcoprime quadraticsで各divisor支持⊂R。
  全βはdiagonal、entriesは`dlog(f/A),dlog(f/B)`。
  conic二型では{A,B}={Q,L²}、三線ではdisjoint-supportの全二degree monomialsを保持する。
* 011：G3=u(a1u+b1v)(a2u+b2v)、a_i,b_iはlinear、各row primitive。
  b1,b2はRのline components、δ=a1b2−a2b1はnonzero quadraticでR支持。
  βの全式は次節で与える。

SFなので各determinantは非零。係数zeroの境界でdegree flag又はpartitionが変わる場合は、
上の別flag又は21/3として保持し、同じopenへ外挿しない。

## 5. SF split011：全collisionと非log部分を保持する式

B=(a1b2+a2b1)/2、C=b1b2と置くと、stable cubicのtwisted GL2 Lie stabilizer零からβは一意であり

```
β11=dlog(f/δ), β22=dlog(f/C), β12=0,
β21=(B dlogδ−dB)/C.
```

全global条件はfβの全係数がquadratic polynomialであること。
この式から次の全場合が得られる。

b1=b2=Lなら、δ=L(a1−a2)=const*L MでMはRのline component。
a=(a1+a2)/2として

`β21=(a dlogM−da)/L`。

M=LではLのgeneric点で二位poleを持ち、primitive rowに反する場合以外消せない。
従ってM≠Lが必要。conic二型はline componentが一つだけなので011は全空。
三線型ではMを他二lineのどちらかに選び、aは全linear formを保持する。
これはtrue nonlog connectionを含む。aを勝手にR-preserving normal form一点へ固定しない。

b1=L,b2=Mが相異なる場合、δがLを含めばa1もLで割れrow1は非primitive、Mも同様。
従って三線R=L M Nでδ=const*N²だけが残る。
δ∈(L,M)なのでtriangleのN(L∩M)≠0では不可能。
concurrentではN=L+Mと尺度を合わせ、common normal translationを使って全解は

`G3=u(Nu+Lv)(−Nu+Mv)`

へ直せる。一般linear解との差はa1→a1+tL,a2→a2+tMだけで、これはactual normal shear。
この全sourceはc（pencil外の第三base座標）非依存であり、β21=dlogL−dlogM。
この段階はfirst connectionの分類であり、actual quadratic P0/higher liftの存在を主張していない。

## 6. Repeated21/3の全primitive degree flags

generic21ではmultiplicitiesが異なるので両rootはC(S)上rationalであり、nontrivial root coverは無い。
全factorizationは

`G3=h ell² m`, `deg ell=d`, `deg m=e`, `deg h=2−2d−e≥0`。

従って全flagsは正確に

`(d,e)=(0,0),(0,1),(0,2),(1,0)`。                      (3)

generic3ではG3=h ell³、3d≤2なので**全triple rootはconstant**でありdeg h=2。
このrationality/degreeの議論をSFのnontrivial coverへ誤用しない。
hとδ=ell∧m（21のみ）のdivisor支持はRに含まれる。
primitive rowのcodimension二basepointsは全て許す。

reduced conic+line R=LQではdegree二支持はQ又はL²、degree一支持はLだけ。
三線R=L1L2L3ではdegree二支持は全Li²及びLiLj、degree一支持は全Li。
marked componentの全選択を、同値を取る前に保持する。

## 7. 21の全character式：double root側だけが固定される

row frameをM=(ell;m)、charactersをθd,θsとすると

`dM+Mβ=diag(θd,θs)M`, `trβ=θd+θs−dlogδ`。

parallel h ell²m/fから

`θd=dlog(f/(hδ))`, `dθs=0`, `θs(E)=e+1`。           (4)

θd(E)=d+1もdegreeと一致する。
両primitive rowsのlocal unit成分を使うと、両charactersはR沿い高々単極。
closednessからそれぞれconstant-residue logarithmic connectionである。
global closed characterは全て

`θ=sum_i ρ_i dlog f_i`, `sum_i deg(f_i)ρ_i=θ(E)`      (5)

で尽くされる。zero-residueの残差はP2上holomorphic一形式となり零だからである。
三root rigidityを21へ流用してθsを根拠なく固定しない。
非log部分は全matrixのoff-diagonal entryに残る。

00ではδはconstant、βはdiagonal。
double側`dlog(f/h)`、simple側は式(5)のEuler一を持つ全charactersである。
hの全finite choicesとsimple characterの全periodsを保持する。

## 8. Moving degree一の両型の全simple-pole条件

01ではell=u、m=b u+L v、hはRのline component。
10ではell=L u+bv、m=v、hはconstant。
LはRのline component、bはLとindependentな全linear form。
frameの明示逆から両型のoff-diagonalは、上下の違いを除いて共通に

`(b η−db)/L`。                                         (6)

01ではη=θs−θd、10ではη=θd−θs。
どちらもη(E)=1、closedである。
Lにηのnonzero residueがあれば(6)はL上二位poleを持つ。
従って**ηのL-residueは零**、すなわちηはR\Lの全closed logarithmic charactersでEuler一。
逆にこの条件ならfβは全quadratic polynomialとなり、flat connectionが出来る。

conic二型ではη=(1/2)dlogQに一意。
三線ではη=ρ dlogM+(1−ρ)dlogN、ρ∈Cの全line。
bは全linear form、constant normal shearでb mod Lの代表を変えられるが、
base Rを保たない任意PGLを使いbを一つの代表へ固定しない。
normalization inverseはconstant normal GL2の逆だけで、bのbasepointがRの交点にあるとはまだ仮定しない。
actual P0のglobal K consistencyが追加条件を課す可能性を残す。

## 9. Moving simple degree二：全B2と全square determinant補集合

02ではell=u、m=A u+B v、A,B quadratics、gcd(A,B)=1、hはconstant。
δ=Bのdivisor支持はR。式(4)からθd=dlog(f/B)。
η=θs−θd、η(E)=2とすると

```
β11=θd, β22=θs−dlogB, β12=0,
β21=(Aη−dA)/B.
```

Bのprime multiplicity m_iでAはunitである。
simple-pole条件はまずres_(f_i)(η)=0を強制する（そうでなければorder m_i+1のpole）。
m_i=2ではさらに(Aη−dA)がそのlineで割れなければならない。
これはtangent derivativeだけでなく、ambient normal derivativeの係数にも課す条件である。

**Bがreduced quadratic。** B=Q又はLiLjで、残りは一line L、f=LB。
ηはLだけに支持されEuler二なのでη=2dlogL。
従ってθs=3dlogL、θd=dlogLであり、全A quadraticを保持する。
`fβ21=2A dL−L dA`は全polynomialである。

**B=L²、conic二型。** ηはQだけに支持されEuler二なのでη=dlogQ。
L-divisibilityはd(A/Q)の全coefficientsがL上零という条件。
tangent方向からA/QはL上constant c、normal方向からA−cQがL²で割れる。
degree二なので

`A=cQ+dL²`, c≠0。

c=0はprimitive性に反する。actual normal shearとscaleでA=Qへ直せる。

**B=Li²、triangle。** 他二lineをM,Nとする。
η=ρ dlogM+(2−ρ)dlogN。
Li上の二交点は異なるため、dlog(A|Li)=η|LiからresiduesはA|Liのzero orders、
従ってρ=0,1,2である。
Q2=M^ρ N^(2−ρ)とすると前と同様

`A=cQ2+dLi²`, c≠0。

全三choicesを保ち、shear/scaleでA=Q2へ直せる。

**B=L²、concurrent。** L=a,M=b,N=a+bに取り
`η=ρ dlogb+(2−ρ)dlog(a+b)`、ρは全C。
第三base座標cを保持して一般quadratic Aを書く。
L上のtangent equationからA|L=α b²、α≠0。
`A=αb²+a(βb+γc)+δa²`とすると、normal derivativeの条件は
`γ=0, β=α(2−ρ)`を与える。
従って全sourceは

`A=α[b²+(2−ρ)ab]+δa²`, α≠0, ρ∈C。

normal shear/scaleでα=1,δ=0へ直せ、全ρを保持する。
この全02 sourceは第三base座標c非依存へ還元される。
ρの特殊値又はAのfactorization境界もこのwhole lineに含めている。

## 10. Triple rootの全closed triangular module

全G3=h u³、hはR支持のquadratic、constant normal GL2でuを選ぶ。
root保存から

`β=[[p,0],[r,q]]`, `2p−q=dlog(f/h)`。

pとqはclosed logarithmic characters、Euler一。
pは式(5)の全Euler一を保持し、q=2p−dlog(f/h)で決まる。
off-diagonal rは**全**H0(Ω1_P2(3))/fの八次元sourceに取り、flatnessは

`dr+(p−q) wedge r=0`。                                (7)

式(7)はfixed periodsでrの全global linear equationである。
全matrixをlogと仮定せず、nonlog rを捨てていない。
これ以上のPGL正常形、actual quadratic P0又はhigher-sourceの次元はまだ求めていない。

## 11. このwhole入口から残る正確な数学

SFの全monodromy coversと全global degree budget、split SFの四flags、
21の全四flags、tripleのconstant rootを全て記録した。
新しい完了は、moving degree一の全simple-pole character条件と、02の全square/reduced determinant分岐を
normal derivativeまで含めて還元したことである。

次の未了事項は、nontrivial SF covers上の低degree sectionsの全分類、
上の全connectionsがglobal quadratic P0を与えるための三normal K係数のconsistency、
full P1三十/P2九と実D既約性への全sourceの橋、及び必要なG上界又は既成分認識。
meromorphic connectionの存在だけでPoisson P0やhigher liftが存在するとしていない。
旧R非被約のfinite profilesをこの被約Rへ数値外挿していない。
新source probe又はGB/jetを始めず、この全数学入口から最弱の次手を選ぶ。
