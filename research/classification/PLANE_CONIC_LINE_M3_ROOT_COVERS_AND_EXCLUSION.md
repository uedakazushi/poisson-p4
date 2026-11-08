# 被約conic+line零曲線・plane位数三：全root coversと横断型の全空性

2026-11-08、PRIMARY G、producer plane_collision_review。
**純幾何producer証明、独立review前。** 新source probe、minor、GB、jet又はCASは実施していない。
actualD被約既約quintic、rank四π、rank二plane S、R=LQ（Q smooth irreducible conic）、normal位数三を対象とする。
全SF root covers、generic21/triple、tangent/transverseを保持する。

**結論。** 横断L+Qの場合は全sourceが空。
接線L+Qの場合は、全sourceのfirst connectionを末表の有限四型へ必要還元する。
接線型のfull higher源のG上界・E帰属は未完了である。

## 1. 既に独立確認されたconstant-root障害と全分岐

`PLANE_M3_CONSTANT_ROOT_OBSTRUCTION.md`、同separate reviewを用いる。
constant normal root uについて、全higher源でu-rowの非可除項はv²C、C constant。
full actualD3のv³係数は

`[D3]v³=C(f)−det(E,C,Y)=3 f q(C)`。                    (1)

root uならこの係数零、C零ならDはu因子を持つのでq(C)=0、C≠0が必要。
R=LQ reduced、Q nondegenerateから全constant-root sourceでq=dlogLが必要である。

既whole degree還元により、21の全moving flags01/02/10は両conic配置で空、
00はh=Qかつβ=dlogL Iへ必要還元された。
triple rootはconstantで、全first connectionsは

* S：β=dlogL I。
* N：β=dlogL I+(dlogQ−2dlogL)N0、N0 nonzero nilpotent。
* D：β=diag((1/2)dlogQ,dlogL)。

へ還元された。N0零はS、nonzeroはconstant normal GL2で一つのnilpotent代表へ移す。
全split SFもSだけに残る。
本稿は、残るnonsplit SFの全C2/C3 coversを処理して、この表を閉じる。

## 2. Universal weight一curl：横断型の全constant-root branches

元全P1三十/P2九から、constant root uの同じbad vector Cに対し

```
Jac_iuv|v³ = k11*C+(DY)C,
Jac_iju|v² = (2Y−X) wedge C−P_{C(f)}.
```

従ってprojective Jacobiのquad(i,j,u,v) v³係数は

`E wedge(k11 C+DY C)−(2Y−X) wedge C+P_{C(f)}=0`。

base volumeでvectorに双対化しcurlを取ると

`(2k11 I+4DY−DX)C=0`。                               (2)

既m二curl証明のliteral identityであり、導出にはk11≠0、R被約又はnormal quadratic SFを使わない。
off-diagonal mixed fieldと他の全higher coefficientsはこのnormal次数へ寄与しない。

conic+lineのclosed horizontal formは一dimで、
`α_V=dlogQ−2dlogL`に対応するtracefree linear fieldをVとする。
S/NではX=Y=V、DではX=−V/2,Y=V。
X,Yは比例しfを保つためD2のuv係数は3f k11で、位数三からk11=0。
(2)はどちらもDV C=0へ帰着する。
横断座標L=c,Q=ab−c²ではV=−a∂a+b∂b、kernelはC∂c。
一方q=dlogLからC(c)=0。よってC=0で矛盾。
**S/N/Dの全横断sourceは空**、Nのoff-diagonal係数もこの議論へ追加のC方向を与えない。
接線座標L=c,Q=ac−b²ではV=2b∂a+c∂bでkernel=C∂aとなり、C(L)=0と両立する。
従って接線型へ空性を外挿しない。

## 3. SF root coverの全型

両conic配置でπ1(P2\(LQ))=Z。
SF degree三root incidenceのimageはS3のcyclic subgroupなので、全coverはsplit、1+2（C2）、connected3（C3）。
S3 connected coverはここには無いが、concurrent三線の別残余では引き続き保持する。
splitは§1で処理済み。

C2のfixed rational rootをprimitive degree dとすると、residual binary quadraticのbase degree e≥1で
`d+e≤2`。
従って全nonsplit場合はd=0又は1。d=2ではresidual coefficientがconstantでsplitとなる。
d=0で全coefficientをdegree二へまとめた場合のcommon factorも、次節で消去又は保持する。

## 4. C2、constant fixed rootの全source

G3=u(Au²+2Buv+Cv²)、A,B,Cはquadratics、C及びT=AC−B²は非零でR支持。
ordinary complement characterはq=dlog(f/C)。
(1)のconstant-root障害よりq=dlogL、従って尺度を合わせC=Q。
Tはdegree四なのでconstant倍のL⁴,L²Q,Q²の全三場合だけ。
nonsplit C2ではTがC(S)上平方でないので、T=τL²Q、τ≠0。
polynomial A=(B²+T)/Qからprime QはBを割り、B=bQ。
constant normal shearでBを零へ戻し、全leading sourceは

`G3=u(τ L²u²+Qv²)`。

全connectionはD型であり、§2から横断型は全空。
接線型ではこの全D profileを残す。common base factorを持つ場合はC=Qと両立せず、別に落としていない。

## 5. C2、moving fixed root degree一の全divisor計算

この場合G3=ell F、ellの係数はlinear、primitive、Fはlinear-coefficient binary quadratic。
disc Fはdegree二でR支持、nonsquareなのでnonzero constant倍のQである。
Fの三coefficientsはindependent（そうでなければdiscのrankは三未満）。
従ってT=(F=0)⊂P2×P1はsmooth double cover of P2 branched along Q、T≅P1×P1。
Hはbase hyperplane pullbackでclass(1,1)、normal P1のclassをF0=(1,0)とする。
div(ell|T)のclassはH+F0=(2,1)。
resultant Croot=F evaluated at ell=0はdegree三でR支持。
全場合はconstant倍のL³又はLQであり、ell divisorのprime supportsはpreimages of L,Qだけ。

Lが横断ならπ^{-1}Lはirreducible(1,1)、ramification curveも(1,1)。
その非負整数和は(2,1)にならない。
よって**moving C2の横断sourceは全空**。

Lがtangentならπ^{-1}L=A+B、classes(1,0),(0,1)、ramification RQ=(1,1)。
class(2,1)となる全effective combinationsは2A+B、又はA+RQだけ。
standard form

`L=c, Q=ac−b², F=a u²+2buv+cv²`

ではA=(u=0) on Tで、2A+Bはell=c*uというnonprimitive rowになり、このdegree flagには無い。
A+RQは唯一のprimitive section、尺度を除き

`ell=b u+c v`, `G3=(b u+cv)(a u²+2buv+cv²)`。          (3)

restriction H0(O_P2×P1(1,1))→H0(O_T(2,1))はinjectiveなので、同じdivisorを持つ別のlinear ellはscalar以外存在しない。
resultantはcQであり、全zero complementsを含めたdivisor classの列挙である。

この唯一tangent connectionをMと呼ぶ。
z=bu+cv,w=uのrational frameのdetは−c。
determinant twistを含むparallel tensorは−G3'/Qとなり、charactersは
θz=dlogc+(1/2)dlogQ、θw=dlogc。
元frameへ戻すと

```
β11=dlogc, β22=(1/2)dlogQ, β12=0,
β21=(b/(2c))dlogQ−db/c.
```

これはglobal simple-pole sourceで、最後のentryはnonlogを含む。
存在するmeromorphic connectionを、未検査のfull higher tensorと同一視しない。

## 6. C3 coverの全compactificationとdegree二pencils

π1=Zのnontrivial C3 characterでは、meridian relationμL+2μQ=0から両branch weightsはmod3で同じ。
root splitting coverのnormal compactificationは

`Σ: w³=LQ ⊂P3_(a:b:c:w)`

である。inverse characterは同じfield coverを与える。
横断型では二A2（局所w³=xy）、接線型では一E6
（平方完成してz²=w³+y⁴）だけを持つnormal cubic surface。
minimal resolution Y→Σはcrepant、H=−KYはnef/big、H²=3、χ(OY)=1。
これらはhypersurface adjunction、Du Val rationalityとP3のline-bundle cohomologyから得る。
terminal三fold分類等の追加分類は使わない。

三root rational mapsを同時resolveすると、全degree budgetはsum Li·H≤6。
deck transformationが三つを巡回するので各degreeは等しく正であり、従って一又は二。
minimal resolution上でroot mapのpencilから固定curveを除いたmovable Cartier classをDとする。
Hodge indexはD²≤(HD)²/3、DはmovableなのでD²≥0。
Cartier parity D²−HDは偶数。
HD=1ならD²=0しか無くparityに反する。
HD=2ならD²=0。
二general membersのintersection零からbasepointも無い。
HD二のpencilがmultiple fibrationならprimitive moving fiberのH-degree一となって不可能。
従って**全root mapsはY上basepoint-free conic-bundle pencils**F、HF=2,F²=0である。
SingΣ上の未解決basepointや非flat normalizationを仮定して捨てていない。

## 7. 全conic pencilはΣ上のlineから来る

E=H−Fと置く。E²=−1、HE=1。
surface RRでχ(E)=1、H(K−E)=−4なのでh²(E)=0、従ってEはeffective。
Σへのpushはdegree一curve、すなわちprojective line ℓである（exceptional componentsのH-degree零も保持）。
Fのgeneral fiberはsmooth connected genus零、h0(F)=2。
Eのsectionの乗法はH0(F)をH0(H)へ単射し、imageはℓを含む全hyperplanesの二次元pencilに等しい。
よって全root mapsは、**Σ上の全linesからのhyperplane pencilの比**に尽くされる。
exceptional fixed partsを零と仮定せず、二次元linear systemの一致から得た結論である。

全linesも直接列挙する。
projection Σ→P2はP3の点(0:0:0:1)を中心とする。この点はΣに無いので、lineはbase line Mへisomorphically写る。
その上でwはlinear、`w³=(LQ)|M`。
M=Lならw=0のramification line ℓ0だけ。
M≠LならL|Mの唯一零点にQ|Mの二零点も重なる必要があり、MはL∩Qの点でのQのtangentである。
接線L型ではそのMはL自身なので、ℓ0以外は無い。
横断型では二intersection pointsの各tangent Mから三liftずつ、計ℓ0と六linesが全てである。
全line族の一般点だけを用いず、この明示cube条件で全linesを尽くした。

## 8. C3の全leading cubics：extra tangentlineを排除する

ramification line ℓ0=(L=w=0)のpencilはw,L。
Norm(uw+vL)を共通base factor Lで割り、全primitive leading cubicはnormal GL2を除いて

`G3=Q u³+L²v³`。                                    (4)

横断座標L=z,Q=xy−z²で、他のlineは例えばy=0,w+z=0。
そのpencil w+z,yのnormをyで割ると

`G3=xz u³+3z²u²v+3zyuv²+y²v³`。

coefficient gcdは一であり、同じroot mapsを持つdegree二coefficient cubicはconstant倍だけ。
discriminantはnonzero constant倍の

`y² z² (xy−z²)²`。

yはRの成分でないので、U上stable parallel tensorに反する。
他五linesもdeck actionとbase R-preserving変換で同じ矛盾を持つ。
従って(4)だけが全C3 leading candidateに残る。

その全ordinary connectionはdiagonalで

`p=(1/3)dlogf, q=(2/3)dlogQ−(1/3)dlogL`。

これはmonomials u³,v³の二weight equations
`2p−q=dlogL, 2q−p=dlog(Q/L)`から一意に得る。
both axesはparallelなので、(1)のcoefficient identityをroot-zero条件なしでも両axisへ適用出来る。
u³ coefficientは、signを除き3f p(C')=C'(f)、C' constant base vectorである。
(4)ではその係数はnonzero constant倍のQでなければならない。
しかしC'(LQ)=kQならmod QでC'(Q)=0、Q nondegenerateよりC'=0、k=0となる。
矛盾。よって**全C3 sourceは両conic配置で空**。
これは有限なroot-cover line列挙の後、fullhigher actualD3 imageを用いた矛盾である。

## 9. 全結論と接線型の残る四first profiles

SFの全coversを§3–8で処理し、21/tripleは§1の全degree還元で尽くした。
横断型では全surviving constant-root profilesを§2で空にしたので、全m三sourceが空である。

接線型L=c,Q=ac−b²では、全sourceは以下の有限four profilesへ必要還元される。

| profile | 全ordinary connection | 全regular parallel G3 source |
|---|---|---|
| S | dlogL I | Q Sym³(u,v) |
| N | dlogL I+(dlogQ−2dlogL)N0 | C*Q u³（uはnilpotent固定root） |
| D | diag((1/2)dlogQ,dlogL) | span(L²u³,Q uv²) |
| M | §5の全nonlog triangular connection | §5式(3)の全scalar line |

N0は一つのnonzero nilpotent代表。N0零はSに含む。
Dの二coefficient零の各補集合も保持し、SFだけをwhole parallel spaceと混同しない。
Mのparallel spaceが一dimであることも確認出来る。
frame z=bu+Lv,w=uでは可能なrational cubic termsは(z³,Qzw²)/L。
元coefficientのpolynomial性をLで比較すると両amplitudesが等しいので(3)のlineだけとなる。

各profileのglobal quadratic normal Kの存在、全P1/P2及びactualDの既約性は別に確認する。
この表はconnectionだけから未知の全higher-sourceを被覆したとの主張ではない。
有限型へ至る全inverseはbase conic/lineのactual PGL、constant normal GL2、pencil basis changesであり、
rational frameは接続とpole証明にだけ使った。全Poisson tensorsをrational gaugeでPGL同一視していない。

この全有限入口が独立reviewされた後、必要なら四whole profilesに限った一回のbounded全30/9検証へ進める。
本稿はそれを実行した記録ではなく、先行の未閉計算を成功へ読み替えていない。
