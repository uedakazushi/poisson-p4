# K4a actual24 親の四つの equivariant coisotropic plane と全section

2026-11-08、PRIMARY E、純手計算。新CAS/point/profileを使わない。
独立受理済み整数literal curve `P(e)` の e=1（全非零parameterの同一orbit）を使う。
全coisotropic planeの分類とはせず、下の一つのregular equivariant branchだけを尽くす。
section divisorの一致から未知full6磁気源の帰属は推論しない。

## 全branchとcoisotropy方程式

K4a P(0) の第二invariant planeは `{a=u=0}`。
torus weights `(c,a,b,u,v)=(2,1,1,0,0)` のもとでこのplaneに正則に戻るgraph chartの
全equivariant planeは

```
a=A*e*c, u=B*e*b+C*e²*c.
```

e=1で `z=a-Ac,w=u-Bb-Cc` とする。
`K4A_SMOOTH_24_PARENT_ACTUAL_CURVE.md` の全整数十括弧から
`{z,w}=Pau-B*Pab+C*Pca-A*Pcu+A*B*Pcb` をplaneへ制限する。
非零となり得る cv,b²,bc,c² の全係数は

```
cv: 6A(A-3B),
b²: 6(B-1)(3B-A),
bc: 36C(B-1),
c²: 18C²+A³(-37+27B)+A²(153-159B)+A(-396+540B)+324-432B.    (1)
```

他のdegree二係数は零。従って(1)が全coisotropy条件である。
A≠3BならA=0,B=1、最後が18(C²-6)なのでC=±sqrt(6)。
A=3BではC(B-1)=0、C=0なら最後は81(B-1)²(3B-2)²。
B=1の場合も最後18C²からC=0。
よって全解は正確に

```
(A,B,C)=(3,1,0), (2,2/3,0), (0,1,sqrt6), (0,1,-sqrt6).     (2)
```

zero/exceptional branchを落とさず、指定点を試した結果ではない。

## Actual anchor determinantとplane sections

baseを(c,b,v)、normalを(z,w)とし、
`U=({c,z},{b,z},{v,z})|plane`, `V=({c,w},{b,w},{v,w})|plane`。
coisotropyによりfullrank4 quinticのplane sectionは、非零定数を除き `det(E_base,U,V)`。
actual tensorからこのdeterminantを手計算すると次の全表になる。

| (A,B,C) | section D の非零定数を除く式 |
|---|---|
| (3,1,0) | `c²*(v-9c)³` |
| (2,2/3,0) | `(3v-19c)*[b²+2c(3v-19c)]²` |
| (0,1,r), r²=6 | `c²*C_r(c,b,v)` |

ここで

```
C_r=3cv²+b²v+2rbcv-12c²v+4rbc²+60c³.
```

例えば最初のplaneは

```
U=(-60bc, -54cv-60b²+486c², -60bv),
V=(-18cv+972c², 18bv+648bc, 36v²+4374c²).
```

Uに60bEを加えると `(0,-54c(v-9c),0)` となり、determinantは
`-2916*c²*(v-9c)³` と直接得られる。
第二planeでは `T=3v-19c` と置くと同じ列操作が
`det=-12*T*(b²+2cT)²` を与える。
最後の二planeでは r²=6 を使って `det=-3888*c²*C_r`。

## 配置と正確な射程

第二planeの Q=b²+2c(3v-19c) はrank3 smooth conic、L=3v-19cはQの接線。
最後のC_rはc=0でb²v、affine c=1でB=b+rと置くと

```
C_r=3v²+B²v-18v+4rB+36.
```

勾配零条件はBv=-2r、B²+6v-18=0。
候補vは `(v-2)²(v+1)=0`、v=-1はC_r≠0。
従って唯一のsingular pointは `[c:b:v]=[1:-2r:2]`、rank一Hessianと非零cubic方向を持つcusp。
c=0はcuspを避け、ordinary tangencyと単純交点を持つ。
これは既nil22 image double_A型の配置であり、残plane1+1+3やstable二配置そのものではない。

固定coisotropic planeを二normal方向へ一様に収縮すれば、実Poisson normalcone familyはpole-free。
ただしその端点のmagnetic residualは制限され、上のdivisor配置だけから全六磁気源へ拡張できない。
plane移動の全defect map又は真の全source同定が別に必要である。
本ノートは四つのactual boundary flagsを保存するだけで、K4a親をT37へ認識していない。
