# moving-simple残源：三fixedのhyperplane因子と全十Ore modules

2026-11-08、producer純証明。旧三十七source/七vertex reportのbytesを保持する。
本noteではGに不要なb-dependent源を展開せず、最弱のnormal-row因子条件を先に使う。
十moduleの計算・全G閉鎖はまだ未了。

## 1. 三fixedはbad normal column零で空

DEG1_HA_CA_OTHERのt=−8/3、DEG1_HB_CA_OTHERのt=4/3、DEG2_CA2_0について、
fullP1必要kernelの末三行、すなわちz-rowのw²C=(w²c0,w²c1,w²c2)の全座標は零。
これは旧37reportのcomplete kernelを全列で見た結果で、特定higher点の評価ではない。
元P0のz-rowは全z倍でKのw²項も零、P2九はbase-baseだけなので、
全actual Πのz-rowがz倍、Euler/Pfaffianからactual Dにz因子がある。
従って三全sourceは既約quinticに矛盾する。

DEG2_CA2_0の真b-dependent Euler障害は数学的に存在し、旧source_covered=falseは変えない。
しかしその障害を解くことはG排除の必要条件ではない。
三fixedの全pointwise Poisson sourceをrank二と呼んだものでもない。

## 2. 残四sourceの全Q/W

残る二whole lineと二fixedはvertex cで全kernel＋FULL P2九をEuler gaugeし、
Π=Q(y)+∂c∧W、y=(a,b,z,w)、全W=W0+十八係数へ入る。
Q-sourceは全u一係数と全central F三係数。
以下のQで書かれ、fullstrictJacobiは恒等零。

HA line：k=−(3t+5)/(3t+8)、

```
Qab=F, Qaz=0, Qaw=u k z², Qbz=u z², Qbw=2u zw, Qzw=0,
W0=(a²,−2ab,−2az,3taw+(3t+5)bz).
```

D3非零はu非零、t≠−5/3,−8/3。二excluded pointは元全leading像零である。

HB line：R=(3t+2)/(3t−4)、S=−(6t+1)/(3t−4)、

```
Qab=u(Raw+bz)+F, Qaz=0, Qaw=u zw,
Qbz=−u Rzw, Qbw=u Sw², Qzw=0,
W0=(a²,−2ab,az,3taw+(3t+2)bz).
```

D3非零はu非零、t≠4/3,−2/3。これも元全leading像零。

TRIPLE fixed：

```
Qab=F, Qaz=0, Qaw=−u z², Qbz=u z², Qbw=u(2zw+z²), Qzw=0,
W0=(0,−3a²,0,−3aw+3(a−b)z).
```

CB2 fixed：

```
Qab=F, Qaz=u z², Qaw=2u zw, Qbz=0, Qbw=−u z², Qzw=0,
W0=(a²,−2ab,4az,10aw−3az−3bw).
```

両fixedもD3非零からu非零。
Q四＋W十八への元全source coordinate mapを一有限jobで確認する。
二lineの旧nonempty exceptional kernelも全十三で、同じcohom＋common-sourceが全列を覆う。

## 3. actual central coboundaryと全zero補集合

a→a+αz+βw,b→b+γz+δwのexact polynomial changeにより、central Fの変化は

```
HA:       u[(kδ−α)z²−2βzw],
HB:       u[−γz²−(S+R)βw²],
TRIPLE:   u[−(δ+α+β)z²−2βzw],
CB2:      u[(γ+β)z²+2δzw].
```

したがってHA/TRIPLE/CB2は全Fをimage span(z²,zw)で割り、F=εw²、ε零/一へ入る。
HBはS+R=(1−3t)/(3t−4)。t≠1/3ならimage span(z²,w²)でF=εzw、ε零/一。
t=1/3ではimageがspan(z²)へ落ちるのでF=εzw+ηw²、ε,η零/一全四を保持する。
η零はwhole lineの値で、η一の二fixedを追加する。

全sourceで実diagonal/tensor scaleによりu一へ戻せる。
HA/HBはπSを保つtensorσ=sa/scとsb=sa sw/szを使い、
HAではu'=σ sb u/sz、HBではu'=σsa u/sz。
残central一係数は自由diagonal characterで零/非零を一へ正規化出来る。
HB t=1/3ではzw/w²の二charactersが独立で、全四zero/nonzeroを覆う。
TRIPLE/CB2ではsa=sb,sw=sz,σ=sa/scがP0を保ち、uと残central係数を同様に戻せる。
逆はminus shearとnonzero diagonal inverseである。
positive-weight shearはP0/leading weight零を保ち、全W十八へ戻る。

## 4. 一回の全十module義務

全残sourceは

- HAのε零/一：二whole t line。
- HBのε零/一：二whole t line。
- HB t=1/3でη一、ε零/一：二fixed。
- TRIPLEのε零/一：二fixed。
- CB2のε零/一：二fixed。

という十whole modulesへ入る。
全full350 affineW十八とactualD全affine係数、全t分母/最大minor/augmented/factor分母例外を検証する。
proper factor又は不整合が最弱の空性。必要なら全original Q-source dimension四とWkernelでG上界を得る。
central zeroではQ-source dimensionを一減出来る。
HA/HB whole fixedplane B=9、fixedt壁B=8、TRIPLE B=8、CB2 B=9。
全projective上界はB+dim(Q-source)+Wkernel+5。

元旧coarse false/kernel/vertex障害は全て保存し、本finite normalizerだけをG閉鎖とは呼ばない。
generic triple、m二、reduced Rは別残余。
