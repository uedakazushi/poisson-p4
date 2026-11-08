# Smooth cubic scroll A25：同じfiber定数の二境界比較と全source正常化

2026-11-08、PRIMARY E、producer純証明。新operator／sourcepoint／CAS計算無し。
**独立review前。** 旧G≤24をE帰属へ読み替えず、全actual sourceの必要regularityを追加で使う。

結論候補は、A25 node全源は空、A25 cusp全actual源は一つのprojective PGL軌道又は空、である。
従って一つのgenuine A25 cusp構造の既知閉parent帰属を得れば、そのparentはwhole A25 sourceを覆う。
現段階で具体成分名、全Poisson分類又は新component数は主張しない。

## 1. 入口と全full-P2 fibers

受理済み `CUBIC_SCROLL_A25_WHOLE_NONDOMINANCE.md`、
`CUBIC_SCROLL_A25_BOUNDARY_SEPARATE_REVIEW.md` と座標／pole erratumを使う。
S=F1、A=2C0+5f。nodeでは両primitive rows、cuspではdenominator一rowだけprimitive。
numerator rowのdivisorial zeroを許す。全R multiplicityと全extensionは保持する。

Y=Bl_S P4=P(F)→B=P2_(a:b:c)でE=2H−B、D'=H+2B。
全incidenceはa*x3−b*x2+c*x0=0、a*x4−b*x3+c*x1=0。
generic c≠0 fiber座標(X,Y1,Z)=(x2,x3,x4)、conic EはQ=XZ−Y1²。
node／cuspの全二rowsはprojective vector fields

`Ti=[[Ai,0],[ri,0]]`

で表し、(Ui,Vi)=(a,b)Ai+c*ri、conic parameter zにli(z)=Ui*z−Vi。
node f=l1/l2、cusp f=l1/l0。
Ωは対応するordered二rowsのdeterminant。πの非零scalarを一回選び、base Poissonをτ=cΩにする。

既G全境界証明が全actual sourceでc∤Ωを示し、L=(c0)上の全binary恒等式から
nodeは22又は31、cuspはdenominator k0のdouble rootだけにした。
node31は既に全排除。node22のfull-P2 zero schemeは2+2、cuspは4という長さ四scheme。
四点が相異なるとかD'のnormalization-flatを仮定していない。

## 2. 一つのflat fiber定数を固定したexact two-form

generic baseではrow matrix

`Lmat=[[−V1,U1],[−V2,U2]]`（node）、
`Lmat=[[−V1,U1],[−V0,U0]]`（cusp）

のSym²でflat projective fiber coordinates(X',Y',Z')を取る。
conic parameterはf=Y'/Z'で、Q'=X'Z'−Y'²=Ω²Q。
R=Q'/Z'²とおく。nodeのD'_fiberはY'=0、cuspはZ'=0。
全connectionの一意性／flat性により、このflat chartで全inverse symplectic formは

```
omega=omega_B+omega_f,
omega_B=da wedge db / Ω(a,b,1)       (c=1 chart),
node: omega_f=(1/k)*dlog(R) wedge dlog(f),
cusp: omega_f=(1/k)*dlog(R) wedge df,
k∈C*.
```

ここでdは**total exterior derivative**。同じ一つのconstant kを全baseで使う。
inverse conventionを全体で反転しても両summandが同時に反転し、以下のcompatibilityは変わらない。
これは旧G proofのgeneric flat P2形を、constantを別divisorごとに吸収せずに固定しただけである。

a=1 chartではomega_B=db wedge dc/(cΩ)、従ってres_L omega_B=−db/Ω_L。
正しいincidence座標ではQのgeneric L-orderは一。
Z'はgeneric p^-1L上unitなのでord_L R=1、fはf_Lへ正則にspecializeする。

c≠0のgeneric reduced Ω上、denominator rowはprimitiveなのでZ'はgeneric fiberでunit、
ord_Ω R=2、fはcollapsed row-ratio f0へspecializeする。
従ってfiber residueはnodeで(2/k)dlog f0、cuspで(2/k)df0。
E,D'の外のこれらのvertical divisorsではomega全体がregularなのでbase residueと相殺する。

## 3. Node22の全actual pairを一つへ戻す

旧全binary／trace還元で両Aiはpure nonzero nilpotent、fixed pointsは相異なる。
actual PGL2とnodeで許された各rowの非零scaleだけを使い

`A1=Nlower, A2=Nupper`

にする。Aut(S)=Aut(P2,p)のblock shear[[I,0],[t,1]]はriをri+tAiへ動かす。
二shearsでr1=(0,R0)、r2=(S0,0)へ戻る。
両primitive rowsによりR0,S0は非零。

block diagonalH=diag(d1,d2,k0)で共役し、両Aiを元へ戻すrow scalesを掛けると
R0→k0*R0*d1/d2²、S0→k0*S0*d2/d1²。
d2=1、d1³=S0/R0、k0=d1²/S0で両方一。
全unitとfinite cubic-root選択、各actual Hの逆が明示され、任意GL2 row mixingを使わない。

よって全node22 actual sourceのpairは

`(U1,V1)=(b,c), (U2,V2)=(c,a), Ω=ab−c²`

へ戻る。以下はこのwhole正常形の必要regularityである。

## 4. Nodeの二境界が要求するkは反対符号

L上、a=1ではcollapsed conic parameter z=b。
f_L=(b*z)/(−1)=−b²、Ω_L=b。
従って

`res_L omega=−db/b+(2/k)*db/b`。

regularityはk=2を要求する。

Ω上はc=1、a=1/bでf0=U1/U2=b。
da wedge db/ΩはdΩ wedge db/(bΩ)だからbase residueはdb/b。
fiber residueは(2/k)db/bであり、regularityはk=−2を要求する。
一つの同じ非零kに矛盾。従ってnode22全sourceも空。
node31の既排除と合わせ、A25 nodeの全zero／collision complementを尽くす。

## 5. Cuspの全denominatorとnumerator

primitive denominatorのAiはpure nilpotent rank一。
適切な2×2 basisではA0=[[0,0],[n,0]]、r0=(R0,S0)、n,S0非零。
したがって全3×3 T0はregular nilpotent、ker T0=p。
columns(e1,T0e1,T0²e1)の行列はdet=n²S0非零でpを固定する。
このactual Aut(S)共役により

`T0=N3=[[0,0,0],[1,0,0],[0,1,0]]`

へ戻る。旧whole cusp還元のfull numeratorはそのflagを保ち、

`T1=[[lambda,0,0],[b0,mu,0],[P,Q,0]], mu*(lambda−mu)≠0`。

許されたpresentation変換ell1→ell1−b0ell0でb0を零にし、P,Qを新係数とする。
第二rowをprimitiveと仮定していない。
このwhole正常形では

```
(U0,V0)=(b,c),
(U1,V1)=(lambda*a+P*c,mu*b+Q*c),
Ω=lambda*ac−mu*b²−Q*bc+P*c².
```

### 5.1 Ωのnonreduced全補集合

c∤Ωは既に全actual sourceで証明済み。
Ωのcomponent Mがmultiplicity m≥2なら、generic c≠0 fiber上でZ'はunit、
R=M^(2m)*unit、fはregular（primitive denominatorが両係数同時にMでvanishしないため）。
omega_f=(1/k)dlogR wedge dfはM上高々単極。
omega_Bは正確にorder mのpoleを持ち、相殺不可能。従ってΩはreduced。
full-P2 fiberとなるcodimension二点を捨てず、generic divisorでの矛盾をpolynomial全sourceへ適用する。

lambda=0ならΩはb,cのbinary quadraticで、reduced二linesとなる。
各lineでrow ratio f0はconstantなのでfiber residueは零、一方base residueは非零。
これも不可能。従ってlambda非零。
mu非零だからΩはsmooth irreducible conic（そのdetはlambda²muの非零constant倍）。

### 5.2 exact constant comparison

L上a=1、f_L=(lambda−mu)/b、Ω_L=−mu*b²。
従ってregularityは

`1/mu−(lambda−mu)/k=0`, すなわち `k=mu*(lambda−mu)`。

Ω上c=1ではa=(mu*b²+Q*b−P)/lambda、row-ratio f0=mu*b+Q。
base residueはdb/lambda、fiber residueは(2mu/k)db。
従って `k=−2lambda*mu`。
同じkを比較して**mu=3lambda**。lambda=0や重複Ωを隠したまま割っていない。

## 6. 残P,Qの全actual逆と一軌道性

numeratorの非零scaleでlambda=1、mu=3とする。
T0=N3、T1=[[1,0,0],[0,3,0],[P,Q,0]]。
actual Aut(S)のH=I+uN3+vN3²で共役すると

```
(T1')21=−2u,
(T1')32=Q+3u,
(T1')31=P+v−3u²−uQ.
```

presentation row-addition T1'→T1'+2uT0で第一のoff-diagonalを消すと、
Q→Q+5u、P→P+v−3u²−uQ。
u=−Q/5、v=−P+3u²+uQで両方零。
Hのdet一、inverseはI−uN3+(u²−v)N3²。
使ったdivisionは非零constantsと既明示lambda,mu,n,S0だけであり、P,Qのzeroを捨てない。

全actual cusp sourceのpresentation pairは正確に

`(T0,T1)=(N3,diag(1,3,0))`

へ戻る。各presentationはker dfとσ、従ってπ_Eをprojective scalarまで一意に決める。
既global vanishingH0(Y,Lambda²TY(−E))=0により元πもprojectively一意。
よって全actual cusp sourceは一actual PGL orbit又は空。
この帰結はseedpointのmembershipからhigherを推測したものではなく、全sourceのactual逆である。

## 7. E帰属に必要な最後の陽性入力

任意の一つのgenuine A25 cusp tensorが既知closed PGL-stable component Kに入ると分かれば、
§6の全逆によりwhole cusp sourceがKへ入り、node空と合わせwhole A25を覆う。
現に別担当genericがtwo concurrent-plane cuspscroll examplesのactual同一軌道とknown-parent mapを検討している。
そのpointwise陽性証明は別receiptで固定し、本純source proofから具体成分名を仮定しない。
旧Gbound／次元だけをE membershipへ読み替えない。
