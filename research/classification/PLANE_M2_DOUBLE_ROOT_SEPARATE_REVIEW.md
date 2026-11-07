# reduced R・m二double rootの接続と全leading還元：別作者査読

2026-11-07、`/root/quartic_review`。**以下の必要還元と表示族を独立に手で確認した。**
対象はrank二Poisson plane、非零reduced zero cubic R、actual D被約既約quintic、
normal位数二、非零double normal leadingである。
全higher sourceの被覆又はG上界は、横断conic+lineの空性以外にはまだ結論しない。
新しいCAS・source matrix・minor・jet計算は行っていない。
全表示式の代数確認は本稿の手証明であり、既証明書の再実行と区別する。

## 1. 全接続入口を保持した査読

`PLANE_M2_DOUBLE_ROOT_CONNECTION_GEOMETRY.md`（1636bb48…）を対象とする。
UFD/二gradingからG2=h(au+bv)²、gcd(a,b)一、a,bの同次degree dは零又は一。
R外のregular flat connectionのparallel tensorは零になれないのでdiv(h)とprimitive rootの共通零点はR上。
非primitiveな境界はcommon factorをhへ戻して扱う。

root equation `d ell+ell beta=theta ell` はprimitive local unitでthetaの単極性を与える。
closednessはgeneric divisorで非log単極を消し、codimension二のroot basepointを理由に余分なdivisorial poleを許さない。
従ってtheta及びtrace betaは全closed logarithmic periodsを保ち、
`2theta−trace beta=dlog(f/h)`、Euler theta=d+1、Euler trace=2。
beta自体の非log entryは除いていない。

moving d一で固定traceの差は `delta beta=omega N_ell`、
N_ell=[[-ab,-b²],[a²,ab]]。gcd(entries)一からf omegaの全係数はpolynomial。
f delta betaのdegree二とN_ellのdegree二から、それらはconstant。
Euler contraction零はconstant-coefficient one-formを零にするので、差は全体で零。

p=(a=b=0)がh上/q外ならh=x,ell=(x,y)と出来る。
tracefree部A dx/x+C dy/xのroot equationをx零に制限すると、
A,Cのresidueは共にstrict upper triangular。
flatnessの二重極係数はC|x0=[A,C]|x0の符号付き式なのでC|x0零。
dy root equationをpで評価してsecond成分一と零が矛盾する。
pはR上であったから、全moving rootでpはq=f/h上にある。

qがideal(a,b)に属すればlinear completion M、det M=qが存在し、
beta=theta I−M^-1dM、theta=(trace beta+dlog q)/2は全simple-pole/global/flat/Euler条件を満たす。
二completionの差はconstant left lower-unipotentなのでM^-1dMは同じ。
これはauxiliary second rowであって、root rowを勝手にmixした同値ではない。
一意性によりこれが全connectionである。

actual P0はPoissonでD(P0)=hQ²、Q=au+bvはrank四、vertex p。
既済 `DOUBLE_QUADRIC_SUPPORT.md` §4の全rank四頂点回避配置の空性をP0へ適用して、pはhにも属す。
この適用は未知のfullhigher Pまでh-Casimirと仮定していない。
結局p in h intersect qの全finite marked配置へ還元される。
constant d零ではbeta lower triangular、omega=(f/h)beta21=N2/h closedを全て保つ。
hの全選択はconic+lineでf/L³、three-linesで111/21/3と全markであり、
未分類のnonlog omega又はperiod零境界をdiscardしない。

## 2. moving d一の全P0正常形とseed帰属

`PLANE_M2_MOVING_DOUBLE_ROOT_NORMAL_FORM.md`（f5c0c5dc…）を確認した。
qが二線ならpのintersection/片側smooth位置を分けてcompletionを上又は下triangularに出来る。
適切なnormal columnがそのnormal variableで割れるため、actualG2の片方のsquare係数は3f kに限られる。
target h L²はf=hLMで割れず、triangle/concurrentの全moving P0は空。

smooth qの横断位置ではh=1、q=1−xz、U=u+zv、V=xu+v。
全scalar periodはX=t(x∂x−z∂z)。全root bracketをLeibnizで戻すと
`{U,V}=qK−txuU−tzvV+3q uv`。
nonzero kに対する{U,V}=kU²はu²係数からq K20=k+txを要求する。
K20はconstantなのでk非零と両立しない。

接触位置だけが残り、q=x−y²、U=yu+v、V=xu+yv、X=t(2y∂x+∂y)。
全normal bracketの多項式性はt非零、k=−t、K=(t+3)u²を強制する。
producer式(2)/(3)の全六括弧を手で再展開した。
A=2t+3、B=t+3と置けば、四strict Jacobiに残る係数は次で全て零である。

* Jxyuのyu係数 −9−3A+6B、v係数 −3B+3t+9。
* Jxyvのu係数は−18(x−y²)+(−3t+3B+9)x−18y²、yv係数9−3A+6t。
* Jxuvのyu²係数2BA−3A−AB+3A−AB、uv係数−6B+9+AB−At。
* Jyuvのu²係数2B²−B²−B²。

Pfは3t(yu+v)²、h-Casimir斉次liftのactual Dは3t h(yu+hv)²。
全t非零を保ち、t=−3等も除外しない。
root frameのmixed bracketsは{x,U}=2tyU、{x,V}=2tyV、{y,U}=tU、{y,V}=tV、
{U,V}=−tU²であり、rational frame上だけの計算も表示polynomial identityとしてq零へ延長する。

seedの既知帰属も最弱に確認出来る。
affine tensorのcomplementary二形式alphaは

`d alpha=2dQ wedge Omega, Q=yu+v, Omega=−t dx wedge du−2t dy wedge dv`、
`Omega²/2=−2t² vol, Pf(alpha)=3t Q²`。
例えばd alphaのdxdydu、dxdudv、dydudv係数はそれぞれ2tu、−2t、4ty。
従って既double-quadric sourceのc=1−k/Aは `1+3/(2t)`。
座標(y,u)をq2の二変数、対応conjugatesを(−2tv,tx)とするとQのlinear radical directionはq2-isotropic。
既受理 `DOUBLE_QUADRIC_SUPPORT.md` §7.3の全isotropic sourceを適用し、
全moving P0(t)はT18/MPS閉成分に属す。
**seedの帰属から未知のhigher sourceの帰属へは推論しない。**

## 3. conic+line：横断全空、接触二必要diagonal型

`PLANE_M2_DOUBLE_ROOT_CONIC_REDUCTION.md`（a476dbae…）を確認した。
constant root uではraw u-columnはu倍、K02零。
conic+lineの全diagonal mixed fields X,Yは一つのtracefree fieldの倍なので、
actualG2のuv係数は3f kappa。targetがsquareなのでkappa零。
fullhigherが既約actualDを持つためのbad P1係数v²CはC非零が必要である。
全off-diagonal sourceを保ったまま、既全projective必要式は

`E cross DY C−(2Y−X) cross C+gradient(C(f))=0`、`(4DY−DX)C=0`。

raw scale `iE i(X−Y)vol=−3f dlog(f/h)` を三成分で確認した。
transverse f=z(xy−z²)ではD=x∂x−y∂y、h=fでX−Y零、h=z³で3D。
tangent f=z(xz−y²)ではN=2y∂x+z∂y、差はそれぞれ零、−3N。

横断h=f、X=Y=lambda Dではlambda非零のcurl kernelはC=c∂zで、
残必要式の第三成分−6czがc零を強制する。lambda零もessential三のfからC零。
横断h=z³ではlambda!=1が同じ排除、lambda=1の全例外式は
`(2bz−cy,2az−cx,2bx+2ay−6cz)=0` なので全(a,b,c)零。
全period/全off-diagonal/全higherでC零になりactualDにu因子を持つ。
movingの既排除と合わせ、**横断conic+line・normal-double全sourceは空**。

接触h=fの必要性はX=Y=−2N、C∈C*∂xだけ。
h=z³ではX=−8N、Y=−5N、C∈C*∂xだけ。
例外lambda=−1は `(cz,−4bz+2cy,cx+4az−4by)=0` が全C零を強制する。
zero periodも別に排除されるので、未処理のresonanceはこの必要還元に無い。
両型のnonlog/off-diagonal/全P1/P2はなお未知として保持する。

## 4. triangle：検出した手誤りの訂正後の全leading還元

査読対象は**訂正版** `PLANE_M2_DOUBLE_ROOT_TRIANGLE_LEADING_REDUCTION.md`
（68859f6e…）だけである。
旧2b1ad43e… §4の「k20非零を一に規格化」は不正確だった。
X=Y=diag(0,−1,1)、Z=diag(a,b,−a−b)では
`det(E,X,Z)=3a xyz`、Z(f)零なので、actualG2=3f(k20−a)u²である。
例えばa一,b零,k20零でもP0はPoissonでactualD非零となる。
producerが本指摘を認めて修正した。旧版にPASS又は実行依存は無く、旧global claimへの影響は無い。

訂正後はk20−a非零を相対normal scaleで一へ戻し、**k20=1+a**とする。
a=−1/k20零境界を全て保つ全二parameter sourceである。
normal shearを存在しない冗長性として使っていない。

全closed off-form moduleも確認した。
h=fはlog二parameter Z。h=x²yはomega=d(Ay/x+Bz/x)+c dlog(x/y)、
h=x³はomega=d(P2/x²)、P2全quadratic modulo Cx²。
x一chartのclosed formからpotentialを取り、horizontal numeratorのdegree二条件が高degree項を消す。
従って高位poleやexact rational boundaryを捨てていない。
h=x³の全Zはx成分零なのでactualG2/u²はideal(y,z)に属し、target gamma x³に反して全P0空。

h=x²yではDelta=diag(−1,2,−1)、Y=diag(s1,s2,s3)、sum si零。
全Pf二係数は独立に

`kappa=s3−s1`、
`G2/u²=3f(k20+B kappa)−xy(Ay+C0x)(s2−s1+4)`。

gamma非零からR=s2−s1+4非零、A零、C0非零、k20=−B kappa。
actual shear v→v−Bu及び相対normal scaleでB零,C0一,k20零へ戻す。
全bad-row三座標の二係数を再展開すると、必要Yは

* C∥∂x：diag(−2s−1,s+1,s)、kappa=3s+1、**s!=−2**。
* C∥∂y：diag(0,2,−2)、kappa=−2。
* C∥∂z：diag(−2,3,−1)、kappa=1。

三Y-lociは互いに交わらず、複数C成分による隠れた枝は無い。
s=−2はleading零となるhigher-order別義務であって、このm二枝で空と判定した点ではない。
従ってtriangleは全moving空/全h=L³空と、二parameter・一parameter・二fixedの四全P0源へ還元される。
fullhigherの成分帰属又はG上界をここから結論しない。

## 5. 境界と保存範囲

本receiptの新結論は、全connection必要入口、moving全P0一parameterとそのseed所属、
transverseconic全higher空性、tangent二必要diagonal型、訂正triangle四全P0源の必要還元である。
constant-root concurrent、tangentの全higher、triangle四源の全higher、nonreducedR、m三等は別残余。
本稿の保存reportは手証明監査を記録し、実行していない新source行列の算術PASSを主張しない。
