# R二重直線＋直線・m三squarefreeの全接続を有限化する

2026-11-07、主義務G。producer純幾何。current code/jobs freeze中に保存する。
全SF normal-root coverのsplit、1+2、connected3を保持する。
全quadratic P0の存在とhigher G上界は、末節の一回のfinite batchで確認する。

## 1. root-coverのc-degreeは零又は二

actual Dは被約既約quintic、S⊂SingDはrank二plane、位数三、
zero cubicはf=a²b、generic normal cubic G3はsquarefreeとする。
chart b=1のt=a/b、r=c/bを使う。
U=Gm×A1のπ1はZなので、degree三root coverのmonodromyはtrivial、C2又はC3。
t=s^d、d=1/2/3のKummer coverで全rootsが分かれる。
coordinate ring C[s,s^-1,r]はUFD、Pic零、unitsはconstant*sの整数power。
primitive root vectorsのpair determinantsとoverall coefficientはunitsである。

各root vectorのr-degreeの和はdeg_r G3≤2。
positive degreeが一つだけなら、他の二つとのunit determinantsによりそのrootもr-independentとなり矛盾。
従って全rootsがr-independentか、degree pattern(0,1,1)である。
C3 coverでは全三degreeが等しいので前者に限る。
C2では後者のdegree零rootがmonodromy固定、splitでもdegree零rootは一意である。

後者のrootはC(t)上へdescendし、primitive homogeneous normal-linear factor ℓ(a,b;z,w)としてglobal G3を割る。
deg_c G3=2なので、もしℓのbase degreeが正なら残るquotientのbase degreeが一以下となりc²項を持てず矛盾。
従ってℓはconstant normal root。actual normal GL2でℓ=zへ置ける。
これで全sourceをc-independent枝とconstant-normal-root/c-dependent枝へ分ける。

## 2. c-independent枝の全一変数接続

stable cubicのtwisted GL2 Lie stabilizer零により、G3/fがfull connectionを一意に決める。
G3がc-independentなら接続もtだけに依存し、full first symbolはEuler gaugeを除き

`P0=∂c∧[a²∂a−2ab∂b+((aM1+bM2)n)·∂n]`, n=(z,w)

となる。ordinary connectionは

`β=(C/t²+D/t)dt, C=M2/3, D=(M1+2I)/3`。

parallel discriminantのcharacter det²からtrC零。
algebraic Riccati rootsのpair differencesのlog derivativeは高々simple poleであり、
`PLANE_DOUBLE_LINE_M3_CYCLIC_CUBIC_FIRST_JETS.md`§3のelementary pole比較でC semisimple非零は不可能。
rootがidentically infinityならRiccati quadratic coefficient零で、残る二finite rootsのlog differenceが
semisimple Cのdouble poleを持つため同じ矛盾となる。

C nonzero nilpotentの場合も同節の比較からRiccati quadratic coefficientは零。
constant normal変換後、finite rootsは

`q=c0/t+d t^e`, c0非零

というlinear Riccati解である。e=−1はnonzero log項を生み不可能。
e零ではlog項の残係数も零でなければalgebraicにならない。
三rootsが全てfiniteならw³係数h=t^p、0≤p≤2で、z³係数のc0³t^(p−3)は
他の指数とe=−1以外ではcancel出来ず矛盾。
従ってinfinityは一つのconstant normal rootである。

残る二finite rootsはsplit又はC2である。
splitではe整数、C2ではeは非整数half-integerで二constantsが互いに反対。
w²z係数のt^pとz³のc0²t^(p−2)からp=2。
他のcoefficientsがpolynomial degree高々二である条件により、splitではe零、C2ではe=−1/2だけが残る。
全actual pairsは、M2をnonzero nilpotent Nへscaleした後

- split：M1=−2I、M2=N。
- C2：M1=diag(−1/2,−2)、M2=N、N=z∂w。

の二つとなる。constant off-diagonal residueはnormal unipotentで除け、同固有値の場合はalgebraic性がその係数零を強制する。
leading amplitudesはfull parallel spaceへ自由に残す。

## 3. M2零の全diagonal pairsは有限集合

M2零ではfinite root monodromyがresidueのprojective nilpotent partを排除するので、
constant normal GL2でM1=diag(λ,μ)へ直せる。
parallel regular monomial

`a^p b^(2−p) z^k w^(3−k)`

の全条件は

`(k−1)λ+(2−k)μ=4−3p`, p∈{0,1,2}, k∈{0,1,2,3}。      (1)

λ=μならλ=4,1,−2の三scalar pairsで、全四normal monomialsを同じbase monomialに掛けられる。
λ≠μの全squarefree spaceは、次のnormal-support pairの少なくとも一つを含む。

`(k1,k2)=(0,2),(0,3),(1,2),(1,3)`。

それ以外の二supportではfixed double rootがあり、squarefree cubicを持てない。
各allowed pairとp1,p2∈{0,1,2}について、(1)の二linear equationsは独立で一つのrational (λ,μ)を与える。
従って高々36 intersectionsと三scalar pairsで全supportを尽くす。
normal swapはλ,μを交換するactual changeなので重複を除ける。
これはparameterの代表点列ではなく、全eigenvalue lineとsquarefree必要supportから得た有限全分類である。
一つのpairに追加monomialsが許される場合、それらも全parallel kernelに保持する。

connected C3の全三pairsは既noteの(3,2),(2,0),(0,−1)として含まれる。
split/C2の全scalar、integer、half-integer resonancesを別に捨てない。

## 4. c-dependent枝：fixed rootと二collision monomials

§1のconstant rootをzへ直し、G3=z(Az²+2Bzw+Cw²)と書く。
A,B,Cはbase quadratic、C及びT=AC−B²はU上nonzeroなので

`C=κ a^e b^(2−e), T=τ a^g b^(4−g)`, κτ非零。

root vectorsのc-degree(0,1,1)とunit determinantからB=B_c(a,b)c+B0(a,b)、B_cはnonzero linearである。
A=(B²+T)/CがpolynomialなのでCはB_c²を割る。
C=abならa,b両方がlinear B_cを割り不可能。
C=a²ならB_cはaのnonzero倍、B0もaで割れるためactual c-shearとscaleでB=acへ直せる。
するとTはa²で割れ、g=2,3,4だけ。
C=b²なら同様にB=bc、g=0,1,2だけとなる。
constant τ amplitudeはconnectionのdlogTから消えるので、全τ非零を含め同じfinite first connectionで扱える。

## 5. full cubic connectionとglobal regularity

G3/fへのparallel conventionを
`d(G3/f)+(βn)∂n(G3/f)−trβ G3/f=0`
とする。fixed root zを保つのでβzw零。
全三remaining coefficientsから一意に

`βzz=dlog f−(1/2)dlog T`、

`βww=dlog f−dlog C`、

`βwz=((B/2)dlog T−dB)/C`

を得る。Euler contractionはI、β−(dlogf/3)Iはhorizontalである。
C=a²,B=acの場合、

`βwz=−dc/a+(g/2−1)c da/a²+((4−g)/2)c db/(ab)`。

g=2,3,4の全てでfβがpolynomial quadratic one-formとなる。
C=b²,B=bcの場合は

`βwz=−dc/b+(g/2)c da/(ab)+(1−g/2)c db/b²`。

fβの最後の係数にa²c/bが残るため、global quadratic Poisson first symbolの必要条件はg=2。
従ってfull candidate connectionsは四つに限られる。

- C=a²、B=ac、T=a²b²/a³b/a⁴の三つ。
- C=b²、B=bc、T=a²b²の一つ。

これらに対しH0Ω1(3)↔sl3の全inverseでmixed columnsを復元し、**full D2=0がraw Kのglobal quadratic解を持つかも確認する**。
meromorphic flat connectionがあることだけからglobal quadratic P0の存在を推論しない。
例えば最後のcandidateではmixed columnsのhorizontal rank二liftがraw Kにbase ratioを要求し得るため、
その全polynomial consistencyを有限証明書へ保存する。
存在する場合はその全P0、存在しなければ明示不整合を保持する。

## 6. finite batchとsource次元の復元

全candidateは§3の有限diagonal集合、§2のnil二pairs、§5の四connectionsで尽くされる。
normalizationにはactual ordered-line base changes、constant normal GL2、tensor scale及びactual c-shearsだけを使った。
Riccatiのrational gaugeは不可能性証明にだけ使い、全Poisson tensorsをPGL同一視するためには使わない。

固定planeのfirst-symbol源の安全なdimension upperは、
diagonal scalarで五、diagonal nonscalarで七、nil scalarで七、nil nonscalarで八、c-dependentで十である。
ordered lines四とπS amplitude一を戻し、normal conjugacy dimensions及びcのactual linear choicesを加えた値。
c-dependentではc mod<a,b>のscaleとnormal nil amplitudeの一冗長性を引き、残るc-choice二を保つ。
stabilizer増大はこのupperを増やさない。

一回のbounded whole-source batchでは全P0のglobal存在/full350/Pf0、
全30/9と全24D3、全parallel-space PDE核、actual leading image/discriminantを保存する。
squarefree leadingが像に無ければ全枝空。
像にある場合はfullhigher dimension上界と上のfirstsymbol/plane/scalarを戻しG≤24を判定する。
足りなければexact whole kernelとoriginal nonlinear Jacobi/actualDを次の構造課題として保持し、
minor/jet/点列を増やすことを進展としない。

R=3Lのm三squarefree全源は別の全二接続/actualD3像証明で空。
本noteの有限化が全計算・独立reviewで閉じれば、nonredR全m三squarefreeのGを処理出来る。
m三repeated root、m二及びreduced reducible Rのm二/三は別残余である。
