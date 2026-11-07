# non-normal CI22：元quadricsの全被約支持分岐

2026-11-07、主義務G、producer plane_collision_review。
全五surfaceのfull150入口と元式は `CLASSIFICATION_NONNORMAL_CI22_TANGENT_20261107.json` に不変保存。
source座標p_iはそのcomplete projective kernel順である。
以下はその全equivalent quadratic rowspaceを使う。元式のschemeを被約同型に置換したとは主張しない。

## 1. regular22_1は二つの五sourceだけを残す

最初の元式(p2−p8)²0からp2=p8=t。
t非零では他の元式からp1=p0=p3=p6=p7=p9=0、p4=−t/2。自由はt,p5の二。
t0ではp1(2p1−p4)=0。
p1非零ならp4=2p1、p3=p6=p7=p9=0、自由p0,p1,p5の三。
p1零なら (p0,p4)*(p7,p9)=0 なので、全supportは次の二linear sourceに入る：

A: p1=p2=p7=p8=p9=0、free(p0,p3,p4,p5,p6)；
B: p0=p1=p2=p4=p8=0、free(p3,p5,p6,p7,p9)。

前二枝は像dimension≤22/23で十分。A/BはAutS三を戻すとcoarse像25を残す。
両全五parameter源のactualDをliteral復元・factorし、既約quintic枝を持つかを別に判定する。
一つの代表点で判定しない。因子が出なければこの補集合は未解決のまま保持する。

## 2. nil32の全sourceはdimension五以下

equivalent元quadricsを保存順E0,…,E16と呼ぶ。
A=p10,P=p4,B=p8,r=p6,t=p3,f=p5,g=p11とする。

A非零：E1,E3,E7からr0,P=3A/2,t=3B/2。
E4は8Af=B²、E8は6p1A+(g+f)B=0。
E9は2A p2+A p9+gB+fB/2=0。
6E10のp0係数は12Aなのでp0を一意に解ける。
従って自由(A,B,g,p2,p7)の五以下。

A0,P非零：E15,E0からr0,t=2B、E4からg=−f。
E2でp1=0、E9でp2=−2p9、E10でp0=0。
自由(P,B,f,p7,p9)の五以下。

A=P0では、E4,E5,E11のcore条件は

`3t²−5r²−3rB−12B²=0, r(t+r−B)=0, 9tB−2r²−18B²=0`。

r0ならt=2B。r非零ならt=B−rかつ(3B+r)(3B+2r)=0。
従って全coreは三line `(t,r,B)=(2s,0,s),(4s,−3s,s),(5s/2,−3s/2,s)`。

s非零の第一lineではE2,E8からg=−f、E6,E12からp1=0,p2=−2p9、E14からp0=0。
自由(s,f,p7,p9)の四以下。
第二lineではE2,E8からf=g0、E12,E6からp2=2p1,p9=−p1/2、
E14からp0=p1²/(2s)。自由(s,p1,p7)の三以下。
第三lineでもf=g0、p2=7p1/3,p9=−4p1/3、p0=10p1²/(9s)。同じく三以下。

s0ではt=r=B=A=P0。
残りは f(g+f)=0、および

`12p1 g+6p1 f−8(g+f)p2−10(g+f)p9=0,`
`6p1²−3p1p9+4(g+f)p7=0`。

g=−f非零ならp1=0、自由(f,p0,p2,p7,p9)の五以下。
f0,g非零ならp2=(6p1−5p9)/4、p7=−(6p1²−3p1p9)/(4g)、
自由(g,p1,p0,p9)の四以下。
f=g0ではp1(2p1−p9)=0、各枝自由四以下。
全反転とその零補集合を列挙したので全source≤5。AutS≥4を戻し像≤24。
各rational mapは必要な全supportを含むsupersourceとして使い、全像がPoissonとは要求しない。

## 3. singularM1の全sourceもdimension五以下

p8非零なら元product式でp0=p1=p2=p3=p7=p9=0、残式でp4=p5p6/p8。
自由(p5,p6,p8)三以下。

p8=0では残条件は

```
p0*p5+p3*p6=0, p1*p2+p6*(p3−2p5)=0,
p1*(p3−3p5)=0, p9*(p3−2p5)=0,
p7*(p0+2p6)=0, p2*(p0+3p6)=0, p1*p7=0, p2*p9=0.
```

p5非零でx=p3/p5とする。p0=−x p6。
x≠2,3ならp1=p9=p6=0、自由(p5,x,p2,p4,p7)五以下。
x2ならp1=0、p2*p6=p2*p9=0。p2零枝は五以下、非零枝は四以下。
x3ならp9=0、p7*p6=p1*p7=0、p1*p2+p5*p6=0。
p7零枝は自由五上の一つ非零式により四以下、非零枝も四以下。

p5零,p3非零ならp6=p1=p9=0、p0*p7=p0*p2=0、最大四以下。
p5=p3=0では p1*p2=0、p7*(p0+2p6)=0、p2*(p0+3p6)=0、p1*p7=0、p2*p9=0。
p2非零ならp1=p9=0,p0=−3p6,p7*p6=0、最大三以下。
p2零,p7零ならfree(p0,p1,p4,p6,p9)五。
p2零,p7非零ならp1=0,p0=−2p6、最大四以下。
全source≤5、AutS≥4により像≤24。

## 4. cone node/cusp

cone_nodeでは全actualD coefficientが元350Jacobi rowspaceに入るため、全Poisson点でD0。

cone_cuspの全二十元quadricsにはp13*p15=0がある。
p13=p15=0なら全actualDは恒等零（full126で確認する）。
p13非零ではp15=0、全元式から

```
p0=p1=p2=p3=p5=p7=p10=p15=0,
p4=2p12, p8=−p6, p14=−p9;
```

free(p6,p9,p11,p12,p13)五以下。
p15非零ではp13=0、

```
p0=p2=p5=p6=p8=p9=p11=p13=0,
p1=2p12/3, p4=4p12/3, p7=(−12p10+9p3)/2;
```

free(p3,p10,p12,p14,p15)五以下。
AutS≥5で全rank四像≤23。零枝も全て保持しており、generic cusp parameterだけの計算ではない。
