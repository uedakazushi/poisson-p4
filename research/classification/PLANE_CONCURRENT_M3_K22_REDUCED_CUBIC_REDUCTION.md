# K22 whole strict-Q源：被約concurrent cubicが許すのは二fixed moduleだけ

2026-11-08、PRIMARY G、producer純証明・独立review前。
新source probe、minor、CAS、GB又はjetは無い。
全二parameter projective row-pairと全Wを保持した係数比較により、
K22枝を二fixed Qとその**whole**W modulesへ還元する。
最後の二moduleのG上界はまだ本稿では証明しない。

## 1. 全K22 sourceとfirst-symbol gauge

先行strict Ore source P=Q+∂c∧W(y)、y=(a,b,u,v)を使う。
K22では

```
Qau=A u², Qav=B u v,
Qbu=C u v, Qbv=D v²,
Qab=−C a v+B b u+C2(u,v), Quv=0,
AD−BC≠0.
```

(A,B),(C,D)は全ordered projective P1 pair、全zero entries/chartを含む。
full W moduleはE_y∧L_WQ=0。
W→W+ell(y)E_yはactual Euler redundancyで、module式を変えない。
ell(a,b)の二coefficientsを使ってW0のbase quadratic部分のdivergenceを除く。
div(ell(a,b)*(a,b))=3ellだから一意に

`W0=(f_b,−f_a,M0(a,b)*(u,v)^t)`

となる。全M0は二2×2 matricesの八係数。
f=f0a³+f1a²b+f2ab²+f3b³は全binary cubic。
reduced concurrent Rはdisc(f)≠0というopen条件であり、normalizerでfを先に一つの代表へ固定しない。

## 2. 元の全moduleから得る四polynomial identities

chart v1でnormal ratioを再びuと書き、p=A−B,q=C−Dとする。
projectivized Q及びc-linear partは

```
pi_ab=−q a+C2(u), pi_au=p u², pi_bu=q u,
Z=(B u,D,0),
P_chart=pi+∂c∧(Y+cZ),
Y_i=W_i−y_i W_v,    y_i=(a,b,u).
```

literal Jacobiから `L_Y pi+Y wedge Z=0`、L_Zpi零である。
符号はminusへ取り違えず、このplusを元Jacobiから用いる。

W0^v=a S_a(u)+b S_b(u)、W0^u=a T_a(u)+b T_b(u)とし、
各S/Tは全linear polynomialとする。

```
R(u)=T_a(u)−u S_a(u), S(u)=T_b(u)−u S_b(u),
Y_u(top degree)=aR+bS,
Y_a(top degree)=f_b−a(aS_a+bS_b),
Y_b(top degree)=−f_a−b(aS_a+bS_b).
```

R,Sは全quadraticである。
full higher WはY_a/Y_bではbase degree≤一、Y_uではbase degree零しか加えない。
C2もbase degree零である。
従って(a,u),(b,u) Jacobiの最高base degree一を比較すると、全solutionに

```
(2p−B)uR−p u²R'−p u²(2f1−2S_a)−q u(2f2−S_b)+qS=0,       (I)
(2p−B)uS−p u²S'−p u²(2f2−S_b)−6q u f3=0,                  (II)
−D R−q uR'+6p u²f0+q u(2f1+S_a)=0,                        (III)
(q−D)S−q uS'+p u²(2f1+S_a)+q u(2f2+2S_b)=0.               (IV)
```

が必要。元fullWをW0だけへrestrictしたJacobiをPoissonと同値と呼ばず、
fullmoduleの次数比較からの必要条件として使う。

## 3. p q非零：reduced fはdiagonal又はanti-diagonal pairへ

(I)のconstant termからS(0)=0。
(II)のu係数からf3=0。
actual swap (a,b,u,v)→(b,a,v,u)でp'=-q,q'=-pだから同じ論法でf0=0。

R=R2u²+R1u+R0、S=S2u²+S1u、
S_a=−R2u+s0、S_b=−S2u+t0と書く。
(IV)のu³係数は−pR2、従ってR2=0。
(II)のu³係数は−A S2。

A非零ならS2=0。
(IV)のu²係数からs0=−2f1、
(III)のu係数からC R1=0、
(I)のu²係数から(p−B)R1=6p f1。
従ってA,C共に非零ならf1=0。
同じactual swapで、B,D共に非零ならf2=0。

いまf0=f3=0なので、f reducedはf1f2≠0と同値。
従って(A=0又はC=0)かつ(B=0又はD=0)。
AD−BC≠0を戻すと全補集合は正確に

* B=C=0、AD非零（diagonal）。
* A=D=0、BC非零（anti-diagonal）。

だけである。p/qのnonzero domains外は次節で全て戻す。

## 4. p0又はq0：radial/axis pointだけが残る

p0ならA=B、det conditionからA=B非零、q非零。
(I)(II)よりS=0及びf3=0。
(IV)からS_b=−f2、(I)からR=−3q f2/Bはconstant。
(III)のconstant termは−DR。

D非零ならf2=0なのでf=a²(f0a+f1b)はnonreduced。
従ってreduced fにはD0が必要で、projective row scalings後

`(A,B)=(1,1), (C,D)=(1,0)`。

さらに(III)からS_a=−2f1となり、W0^vのu係数は全て零。
すなわちglobal W0^vはvで割れる。
q0の場合はactual swapの同じradial/axis pointである。
p=q0ならrowsがproportionalでdet零なのでsource外。

## 5. diagonal/radial-axisはfull higher込みでactual factorを持つ

### diagonal

B=C0、A,D非零。
§3の(I)constant及び(III)constantからS(0)=R(0)=0。
global W0^uはuで割れる。
full quadratic W^uのうちu非可除な残りは、normalweight条件からλv²だけ。
chart v1でu0ではY_u=λ、baseにconstantである。

このchartのpiは

`pi_ab=D a+C2(u), pi_au=A u², pi_bu=−D u`, `Z_b=D`。

元fullmoduleの(b,u)係数をu0で取ると左辺L_Ypiは−Dλ、
Y wedge Zは−Dλなので−2Dλ=0。
従ってλ0、full W^u及び全Qのu-columnはuで割れる。
full Pのu-columnもu倍、actualDにu因子を持つ。
C2全三coefficientsと全higherを保持した結論である。

### radial/axis

p0、(A,B)=(1,1),(C,D)=(1,0)では§4からW0^vはvで割れる。
残りの非可除項はλu²だけ。
normal u1、base swapしたchartを用いると

`pi_a'b'=C2(v), pi_a'v=−v², pi_b'v=0`, `Z_b'=1`。

v0でY_v=λconstant。
fullmoduleの(b',v)係数はL_Ypi零、Y wedge Z=−λなのでλ0。
従ってfull v-columnはv倍でactualDはvで割れる。
q0の対称radial/axisも同じ。
これら全源は被約既約quinticを持たない。

## 6. 全残源は二fixed anti-diagonal Qとwhole W

唯一残るrow pairはA=D0、BC非零。
row scalesでB=C1へactualに戻す。
全Qは

```
Qau=0, Qav=uv,
Qbu=uv, Qbv=0,
Qab=−a v+b u+C2,
Quv=0.
```

全C2 shift imageはspan(u²,v²)なので、actual a/b-normal shiftsでそれを除く。
residual γuvはtrue common tensor/base scaleでγ0/1。
ゆえにwhole残源は正確に二canonical Qε

`Qε_ab=−a v+b u+εuv`, `Qε_av=Qε_bu=uv`, ε=0,1

と、その**全**quadratic W source、module E_y∧L_WQε0へ入る。
元rank二plane/normalorder三/f reduced条件は全Wへ戻す。
Wを選んだpointで置換せず、全first四＋normal八＋higherを保持する。
この二moduleのactualD因子又はambient像上界はまだ未計算・未証明。
次の一回のfixed linear module計算、又はwhole hand module proofが必要である。
本稿だけでK22 whole G/Eを閉じたとは呼ばない。
