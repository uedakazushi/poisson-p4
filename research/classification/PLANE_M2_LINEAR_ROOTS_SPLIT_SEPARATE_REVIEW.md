# m2 linear root rows・split determinant：一Pf係数による全空性の独立検証

2026-11-07、独立 reviewer `/root/quartic_review`。**PASS（純証明）**。
producer `PLANE_M2_LINEAR_ROOTS_SPLIT_EXCLUSION.md`、SHA256
`8efb500a30c7836c8a615afedac1081bea3e6092f4cbeaedee7b85898cced627`
を全文読解し、全正常形・全U-family・actual Pf係数の含意を独立導出した。
新しいsource、線形代数、Jacobi、minor、GB、jet計算は無い。

対象は reduced非零R、restriction rank2 plane、m2、squarefree normal quadratic、
actualD被約既約、primitive linear root rows のflag(1,1)、detMが相異二直線に分裂する全枝。
**この全sourceは空**。conic determinantの全MPS/空性は別の独立receiptで受理済みである。
nonreducedR、他flags、double normal quadratic、m3へ適用しない。

## 全parameter上でscalarなroot-line connection

既独立受理のroot-line class lemmaによってη=θ1−θ2はhだけにpoleを持つclosed log形式、
Euler contractionはd1−d2=0。hは一つのlineなのでη=0。
従って全connectionは

```
β=θId−M^-1dM.
```

θの全residue parameterを固定していない。
M^-1dMのnonlog entryも全て保持する。
このscalar性はreducedRとprimitive rowsに依存し、nonreducedRへ流用しない。

## β_vu=0から全sourceの一Pf係数へ

全normal-order-zero源はbase Jacobian代表、全四mixed linear fields、
全normal-normal quadratic K=k20u²+k11uv+k02v²を持つ。
projective Euler correctionのnormal-linear四成分を使って各mixed fieldをtracefreeにできる。
この補正はbase partを変えず、Kの三係数内で吸収されるので全sourceを保つ真のgaugeである。
tracefree Vの写像 `V→iE iV vol` はinjective。
よってβ_vu=0は対応mixed field V_vu=0を強制し、raw列P_(i,v)は全てv倍となる。

actual invariant `E∧Π0²` をv=0へ制限する。
mixed-v列が零、Ev=0なので、base-base PfとK以外の項は全て消える。
normal u²係数は

```
(Ebase f) k20 = 3 f k20
```

（base Jacobianを1/3で正規化する慣例では f k20）だけ。
必要なのは**非零定数倍のfで必ず割れる**という不変な条件である。
k20はglobal quadratic tensorのnormal-normal係数なのでbaseに依存しない定数。
Euler gaugeをstabilizerと数えたり、任意raw代表にdiv0を課していない。
全P1/P2を戻す前にこのnormal-order-two必要条件が既に破れれば、元の全sourceも空となる。

## image rank三の全正常形

W=im(M)がMat2内三次元ならannihilator Kはrank1。
rank2ならdet|Wはsmooth conicとなりsplit仮定に反する。
rank1 Kのrow relationで一方のrow coefficientが零なら、そのroot rowが一つのconstant normal directionに入り、
degree一の共通因子を持つためprimitiveでない。
従って両係数非零、legitimate row rescalingとnormal右GL2だけで

```
M=[[x,y],[x,z]],  c=z−y, detM=xc
```

に戻せる。marked rootsをmixするleft GL2をsource同値として使わない。
x,y,zは独立linear coordinates、全非零尺度とhの全positionは残る。
直接微分すると

```
M^-1dM=[ [dx/x, (zdy−ydz)/(xc)], [0,dc/c] ].
```

よってβ_vu=0。
target G2=κh(xu+yv)(xu+zv)、κ≠0のu²係数はκh x²。
f=hxcでcとxは独立、hはreducedR条件でx,cと異なるので、これはfで割れない。
上のwhole Pf係数必要条件と矛盾し、全sourceが空。
θとupper off-diagonalのnonlog項を捨てずに結論する。

## image rank二の全U-family

W二次元、detの二root相異なるregular matrix pencilは

```
M=U diag(x,y) V,  U,V∈GL2(C)
```

となる。invertible memberを一つ選び、そのinverseと他の元の積の相異固有値で対角化すればよい。
baseの二linear forms x,yも独立。normal右GL2でVだけを除き、**Uは全parameterのまま残す**。
各rowprimitiveなのでUの四entryは全て非零。
一entry零なら残る同rowをx又はyが共通に割り、別のdegree flagへ戻るためである。

ΘscalarなのでUを除かなくても `β=θId−diag(dlogx,dlogy)`。
特にβ_vu=0、actualu²係数はf倍。
一方targetはκh U11U21 x²、f=h det(U)xy。
yはxと独立、κU11U21≠0なのでfでは割れず、全U-familyが空。
special Uをsamplingしたり、left GL2でparameterを消していない。

image rank一では全linear entriesが同じlinear factorを持ちprimitive条件に反し、
detはsquareになるので今回のreduced split仮定にも入らない。
以上で全possible image ranksを尽くし、codimension二のrow-basepointも除去していない。

## 合成と残余

別receipt `PLANE_M2_LINEAR_ROOTS_CONIC_SEPARATE_REVIEW.md` と合わせ、
**reducedR、m2、squarefree normal quadratic、flag(1,1)の全sourceはempty又はMPS40**。
conic determinantは両rowbasepointsを通るhだけ全SNC曲線、他positionは空。
split determinantは今回の全空性である。
flags(0,1)/(0,2)、nontriangle(0,0)、normal double quadratic、nonreducedR、m3は保持する。
全rank2planeや全E/Gの完了をこの一flagの結果と混同しない。
