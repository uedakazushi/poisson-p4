# Stable three-lines全六磁気源：T03真period曲線による全点被覆

2026-11-08、PRIMARY E_multiple_hyperplane。新しい純actual全source曲線。
本稿作成時には新operator／point／rank／minor／profile／jet／CASを実行していない。
後続の同じ全変数恒等式の有限readback／独立reviewを区別する。
**以下の全恒等式が確認されればstable_three_lines全六磁気源は既知閉成分T03に入る。**
line_conic又は全nil22image laneまでは主張しない。新成分増分は零。

## 1. 真のT03四period・一edge源と全Euler gauge

公開exact sourceは `reports/MPS_S61_RAW_CATALOGUE_20261103_v62/T03.json`。
一edgeの全parameter sourceそのもので、source座標を(X,Y,H,Z,W)=(x0,x1,x2,x3,x4)とする。
別のT03 completionファイルを仮定せず、公開全括弧と全projective Jacobi証明を使う。
H-Casimir gauge後の全六括弧は

```
PXY=aXY+tZW, PXZ=bXZ, PXW=cXW,
PYZ=-(b+d)YZ, PYW=(d-c)YW, PZW=dZW.          (1)
```

元全四periodとのactual可逆線形対応は

```
a=p3, b=-p0, c=p0-p1/2+p2/2+p3, d=2p0+p2+p3;
p0=-b, p1=d+a-2c, p2=d+2b-a, p3=a; t=t0.
```

raw sourceで `Xraw=Praw#dH/H=(0,aY,0,bZ,cW)` とすれば
`P_H=Praw+E_source wedge Xraw` の全H行は零、残りは厳密に(1)。
最高部だけの同一視ではなく、全period／edgeを保持した全十括弧のraw恒等式である。

actual affine Pfaffianは

```
Pf(P_H)=dZW*F, F=(a-b-c)XY+tZW.              (2)
```

FはPXYの式そのものではない。cross termの積を保って(2)を得る。
全sourceのPGL飽和閉包がnamedT03であることは既受理MPS source theoremを用いる。

## 2. 三つの真period derivativeとpole-free実曲線

A,B,Cを自由とし、source period/edgeを**同じepsilon上で**

```
d=1, a=epsilon*A, b=epsilon*B, c=1+epsilon*C,
t=1+epsilon*(B+C-A)
```

と選ぶ。k=a-b-c=-1+epsilon*(A-B-C)なのでt=-kが全epsilonで成立する。
quadricは `F=k(XY-ZW)`、d*k!=0の稠密curve openで相異なる三Hとrank四既約quadricを持つ。
kの定数は-1なので、どのA,B,Cでもこのgeneric openは非空である。
periodとsmoothingの積関係を独立係数へ緩めていない。

sourceからtargetへの実変換を

```
H=h, Z=z, W=w,
X=h-w+epsilon*q1, Y=h-z+epsilon*q2,
tensor scale=epsilon
```

というinverse式で指定する。
forwardは `q1=(X-H+W)/epsilon,q2=(Y-H+Z)/epsilon`、h=H,z=Z,w=W。
source順(X,Y,H,Z,W)、target順(h,q1,q2,z,w)のdetはepsilon^-2。
epsilon!=0でactual PGLと非零tensor尺度である。

X,Yを上のinverse式として、全pushの非零括弧は**正確に**

```
{q1,q2}=A*X*Y+B*X*z+C*Y*w+(B+C-A)*z*w,
{q1,z}=-z*w+epsilon*B*X*z,
{q1,w}=(1+epsilon*C)*X*w,
{q2,z}=-(1+epsilon*B)*Y*z,
{q2,w}=z*w-epsilon*C*Y*w,
{z,w}=epsilon*z*w.                           (3)
```

h行は全零。全空間次数は二、epsilon次数は高々二である。
全normal-pair poleの相殺は、n1=X-H+W,n2=Y-H+Zに対し

```
{n1,n2}=aXY+bXZ+(c-d)YW+(t-d)ZW
       =epsilon*[A*X*Y+B*X*Z+C*Y*W+(B+C-A)*ZW]
```

からliteralに従う。generic source pointでの零からの外挿ではない。
(3)は全epsilonで多項式だから、元sourceの全Jacobiはepsilon零へも延長する。

端点は

```
M=[[-zw, w(h-w)],[-z(h-z), zw]],
R_end=A*h²+(B-A)*h*z+(C-A)*h*w.
```

## 3. 定数actual frameで元stable anchorへ戻す

次のdet729の定数GL5を使う。

```
(Q1,Q2,x,y)=(3q1,3q2,9w-3h,3h-9z), h=h.
```

inverseはq_i=Q_i/3、z=(3h-y)/9、w=(x+3h)/9。
全四mixed bracketsは正確に

```
{Q1,x}=-(x-6h)(x+3h)/3,
{Q1,y}={Q2,x}=-(x+3h)(y-3h)/3,
{Q2,y}=-(y-3h)(y+6h)/3,
{x,y}=0.
```

これはaccepted stable_three_linesの元anchorである。
magnetic端点は

```
R0=r0*h²+r1*h*x+r2*h*y,
r0=3(A+B+C), r1=C-A, r2=A-B,
A=r0/9-r1/3+r2/3,
B=r0/9-r1/3-2r2/3,
C=r0/9+2r1/3+r2/3.                           (4)
```

full frameのdetは729/epsilon²、raw sourceとの違いは§1のEuler項を同frameで移したものだけ。
追加Euler近似又はnonlinear coordinate変更は無い。
source signed Dの変換則はtensor epsilon²とGL determinantにより全体factor729となる。
epsilon零の実Dは厳密に

```
-h²*(x+3h)*(y-3h)*(x-y-3h),
```

元stable fullsourceと一致する。

## 4. 全六磁気をparameter無反転で戻す

targetの任意全六magnetic sourceを

```
R=m0h²+m1hx+m2hy+m3x²+m4xy+m5y²
```

とする。旧全polynomial normalizerの三invariantは

```
I3=m3+m0/27+m1/9+2m2/9,
I4=m4-m0/27+2m1/9-2m2/9,
I5=m5+m0/27-2m1/9-m2/9.
```

次の全period derivativesを選ぶ。

```
A=I3-I4+I5,
B=I3+2I4+4I5,
C=4I3+2I4+I5.                                (5)
```

対応する(4)の係数は

```
r0=m0+18m3+9m4+18m5,
r1=m1+3m3+3m4,
r2=m2-3m4-3m5.
```

最後に同じbase上の有限actual translation

```
Q1'=Q1+3m4*x+3m5*y,
Q2'=Q2-3m3*x                                 (6)
```

を行う。逆は同じshiftを引く。
magnetic変化は全式で

```
-3m3*M11-3m4*M21-3m5*M22
= m3x²+m4xy+m5y²-(3m3+3m4)hx
  +(3m4+3m5)hy-(18m3+9m4+18m5)h².
```

従って(6)はR0から任意のRへ厳密に移す。
(5)と(6)は全六mの多項式であり、一つの磁気parameterも反転しない。
map `m -> (A,B,C,3m4,3m5,-3m3)` のdetは-3。

全inverseも明示できる。shift係数を(s,t,u)=(3m4,3m5,-3m3)とし、

```
m3=-u/3, m4=s/3, m5=t/3,
A'=A-m3+m4-m5,
B'=B-m3-2m4-4m5,
C'=C-4m3-2m4-m5,
m2=A'-B', m1=C'-A', m0=3(A'+B'+C').
```

全zero residue、I全零、mの零係数境界を含む真の可逆六parameter移送である。

## 5. Whole membershipの最弱論理

任意に固定した全六mに対し、(5)でA,B,Cを選び、(3)、定数frame、(6)を合成する。
epsilon!=0かつk!=0の非空稠密curve openではtrue T03 sourceのactualPGL像なのでnamed閉成分T03内。
全六括弧はepsilon零まで多項式で、端点は元full6magnetic stable_three_linesと一致する。
端点はD非零でrank四なのでprojective tensorも非零。
従ってT03の閉性から、**全六磁気源の全点がT03**。
PGL不変性から全飽和とその閉包もT03に入る。

単一磁気代表の所属を全源へ推測したのでも、known boundaryからmembershipを逆推論したのでもない。
local smoothness／unique component／normal支持／newcountは要求しない。
残るstable_line_conicは本sourceと配置が異なるので、この証明へ含めない。
本稿の全source恒等式は次のONE有限jobで確認し、別作者review後に正式採用する。
