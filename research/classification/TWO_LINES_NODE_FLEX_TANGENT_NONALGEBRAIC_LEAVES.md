# node-flex/tangent-through-flex：全磁気商の一般葉に対するresidue障害

2026-11-08、PRIMARY E_reduced_nonnormal、担当 `/root/exceptional_coverage`。
**純手計算。新CAS、source点、operator、rank、profile又は変形jetは無い。**
T37側のbounded invariant-surface条件とproper specializationは別担当の証明である。
本稿はtarget側の全磁気parameterを保つ有理正規化と一般葉を扱う。

## 1. 全六磁気の正則正常化を再使用する

accepted全源と独立査読済みの有限shearを使う。座標(h,x,y,z,w)、h-row零、
既約cubic C=h y²−x(x−h)²に対して

```
U=2x(x−h)∂x+y(3x−h)∂y,
V=−2xy∂x+[2x(h−x)−3y²]∂y,
d=x−h,
P=U wedge ∂z+V wedge ∂w+(−d*w+Ih²+Jhy)∂z wedge ∂w.
```

ここでI,Jは任意の複素数である。正規化前の全R6について
I=r0+r1+r3、J=r2。元公開integer sourceからのconstant frame後ではrは元の四倍である。

全六係数からのfinite actual shearは

```
ell=(r4/2)x+(r5/3)y,
k=−(r1+r3)h+(r3−2r5/3)x,
z'=z+ell, w'=w+k.
```

有限push R'=R+Vell+(d−U)kはR'=Ih²+Jhyを与え、逆はshiftの符号反転、det一。
source係数で割らず全R6を保つ。J≠0へのparameter化lambda=I/Jは本稿には不要で、
J=0,I≠0とI=0,J≠0も同時に扱う。

## 2. cubicの正規化と全restriction

h=1でCはy²=x(x−1)²。
正規化parameterをtauとし

```
x=tau², y=tau(tau²−1)
```

と置く。nodeはtau=±1、flex側の無限遠はtau=∞。
U,Vをこの正規化へ引き上げると、f=tau(tau²−1)に対し

```
{tau,z}=f,
{tau,w}=−tau*f,
{z,w}=−(tau²−1)w+I+J*tau(tau²−1).
```

第一二式はdx/dtau=2tauと{x,z}=2x(x−1)、{x,w}=−2xyから得て、
全式は多項式なのでtau=0での見掛けの除算も恒等的に戻る。
第三式は元の全磁気Bを代入したもので、I又はJを落としていない。

tau≠0,±1の開ではrankは二である。
順序(tau,z,w)の核一形式をfで割ると

```
eta=dw+tau*dz+[J−w/tau+I/(tau(tau²−1))]dtau.
```

さらにtauで割り、q=z+(w+I)/tauと置けば

```
omega=eta/tau
     =dq+[J/tau+I/(tau²−1)]dtau
     =dq+J*dlog(tau)+(I/2)*dlog((tau−1)/(tau+1)).       (1)
```

I/tauのrational exact項をqへ含めることが必要である。
この変更は葉の代数性を調べるためのbirational座標であり、actual PGL変形と称さない。

## 3. 代数的葉があればfinite cover上でexact/logのresidueが矛盾する

s=wとすると新三座標(tau,q,s)の全括弧は

```
{tau,q}=0,
{tau,s}=−tau*f,
{q,s}=I*tau²+J*f.
```

特にHamiltonian vector field of tauは非零単元倍の∂sである。
一般tauを支配するirreducible algebraic invariant surface Wがあると仮定する。
C(tau)へ局所化するとheight一primeはirreducible F(q,s)により生成される。
∂s-invarianceはFが∂sFを割ることを意味し、s次数を比較すると∂sF=0。
従ってFはqだけの非定数多項式で、qはC(tau)上代数的となる。

対応する代数曲線の滑らかな射影正規化Gammaを取る。
tau:Gamma→P1は非定数ゆえfiniteで、qもGamma上のmeromorphic functionである。
葉の接方向で(1)が零だから、Gamma上で

```
dq=−J*dlog(tau)−(I/2)*dlog((tau−1)/(tau+1)).         (2)
```

左辺はmeromorphic functionのexact differentialなので、全点でresidue零。
右辺のbase P1上のresidueは

```
tau=0: −J, tau=1: −I/2, tau=−1: I/2, tau=∞: J.
```

finite cover上でresidueは正のramification index倍となり、標数零では消えない。
I又はJが非零なら少なくとも一つで矛盾する。
従って

```
(I,J)≠(0,0)なら、cubic Cの一般点を通る代数的invariant surfaceは存在しない。
```

非整数／無理数等のgeneric residue条件は使っていない。
I,Jが整数比又は有理比であっても、logではなくmeromorphic qのexact differentialと比較するので同じ障害が成立する。

## 4. 全補集合とT37への最弱適用

I=J=0ではomega=dqとなり、有理first integral qを持つ。
この閉磁気補集合を非代数的と称したり、親から排除したりしない。
本稿はそのprojective葉の次数又は帰属を計算していない。

tau=0,±1,∞の個別fiberも全supportから削除していない。
それらの上にはspecial invariant surfaceがあり得るが、一般cubic点の被覆には使えない。
本稿の開は葉のgeneric algebraicityを検査するためだけに使い、
全Poisson源をその開だけへ置き換えたわけではない。

D=h*x*Cの一般C点はh*x≠0の開にある。
そこを通る、Dに含まれるinvariant surface componentはCに含まれねばならない。
従って(I,J)≠(0,0)では、どの有限次数のeffective invariant二cycleも一般C点を通らない。
別component又はmultiplicityを加えてもこの欠如を補えない。

別担当がwhole T37のrank四境界にdegree≤5 invariant surface cycleによるDの被覆を証明すれば、
この全(I,J)非零源はT37の外となる。元全R6及びPGL像へは実normalizerで戻せる。
閉包のI=J=0境界まで陰性を延長しない。
また全sourceの稠密開がT37の外であることは、その全sourceがT37に含まれる可能性を排除するが、
他のknown親、新成分性又は全E完了は別の義務である。
