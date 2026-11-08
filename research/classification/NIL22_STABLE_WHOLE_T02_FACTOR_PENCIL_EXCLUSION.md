# Stable二全源はT02のwhole factor-pencil閉包の外

2026-11-08、PRIMARY E_multiple_hyperplane。純手証明、独立読解待ち。
新operator、point、rank、minor、profile、normal jet又はCASは実行しない。
全T02 sourceの四periodとedgeを保持する。T03、T16/T17、T36/T37への同じ陰性を主張しない。

## 1. true T02からの全source表示

`MPS_S61_RAW_CATALOGUE_20261103_v62/T02.json` の全sourceを
`(X,Y,Z,W,H)=(x0,x1,x2,x3,x4)` と読む。
H-Casimir Euler gauge後の全六括弧は

```
PXY=F=aXY+tH², PZW=dZW,
PXZ=bXZ, PXW=cXW, PYZ=-bYZ, PYW=-cYW.       (1)
```

元periodとの全域可逆な線形対応は

```
a=p1, b=(p1-p3)/2, c=(p1+p2)/2, d=p0+(p2+p3)/2, t=t0;
p1=a, p3=a-2b, p2=2c-a, p0=d-c+b.
```

raw Pに対し `Xraw=P#dH/H` は
`(0,p1Y,(p1-p3)Z/2,(p1+p2)W/2,0)`。
`P_H=P+E wedge Xraw` を各八entryへ引き算して(1)を得る。
従ってperiodを削ったsubfamilyでもdivisorだけの推定でもない。

一般Dは `H*Z*W*F`、quadric Fはrank三である。
四period全体でPoisson、PGL飽和閉包がnamedT02であることは既受理completed-source theoremを使う。

## 2. 二つのclosed logarithmic-pencil条件

Poisson factor flag `(H,Z,W,F)` を付け、degreeを(1,1,1,2)、D=HZW Fとする。
H-Casimir代表で全体に正則なlinear logarithmic fieldsを

```
Y_Z=P_H#dZ/Z, Y_W=P_H#dW/W, Y_F=P_H#dF/F
```

と定義する。最後のものはprojectively `π#dlog(F/H²)` である。
source (1)では U=X∂X-Y∂Yとして正確に

```
Y_F=aU, Y_Z=bU-dW∂W, Y_W=cU+dZ∂Z.
```

従って一般sourceは次の二条件を満たす。

```
(CZ) some [lambda:mu] in P1 makes lambda*Y_Z-mu*Y_F
     projectively zero on W;
(CW) some [lambda:mu] in P1 makes lambda*Y_W-mu*Y_F
     projectively zero on Z.
```

実際 `[lambda:mu]=[a:b]` 又は `[a:c]` で、場はそれぞれ
`-ad W∂W`、`ad Z∂Z` となる。
同時にa,b,cが零になる境界で比を無理に評価せず、一般flagのclosureをP1係数も含めて取る。

H-Casimir liftはHの各係数chartで一意かつ正則である。
Z又はFによる正確な多項式除算は、非零formによる乗法mapが一定rankのsubbundleを与えることから正則。
これはZ/Hの値又はquadric discriminantを反転する操作ではない。
Fが可約／平方になりHを含む場合も、P_H#dH=0と各因子のPoisson接性から分子のF割切りは残る。
従って全reduced flag incidence上でY_Z/Y_W/Y_Fは正則に延長する。
projectively zeroはE wedge fieldを指定hyperplaneへ制限した全係数零なので閉条件。
P1係数とfactor projective空間はproperであり、(CZ),(CW)の存在条件は全rank四factor-flag closureに保持される。

## 3. 一般hyperplane-collision補題

別の必要条件も全源で用意する。
正則rank四quadratic projective Poisson族で、相異二marked Poisson超平面が同じhへ合体するとする。
一つをregular PGLでhに固定し、他方を

```
H2=h+beta(epsilon)*ell_epsilon, beta(0)=0,
```

とする。betaは差の全係数の最小valuation、ellはh成分零でell0非零。
h=1のquadratic affine代表P_epsilonは正則である。
H2のPoisson接性から各affine座標xjに対し

```
{ell_epsilon,xj}= (1+beta*ell_epsilon)*B_j,
```

B_jはaffine-linearとなる。
B_jはgeneric fieldでの商だが、constant/linear係数を順に比較すればregularである。
具体的にはBの定数は左辺の定数、線形係数は左辺線形係数からbeta*ell*Bconstを引いたもの。
従ってepsilon零で `{ell0,xj}` はaffine-linear。
最高二次bivector Qtopについて `Qtop#dell0=0` が必要となる。
有限base changeや任意の衝突orderでも同じであり、新Taylor orderを選んだ条件ではない。

両stableのoldh=1最高部は

```
Mtop=-p*p^t/3, p=(x,y), Pxy=0,
Pq1q2=Rquad(p).
```

ell=A·q+B·pがQtopのlinear Casimirなら、base成分は
`(A·p)*p/3` なのでA=0、続くfiber成分 `-p*(B·p)/3` からB=0。
従って**全磁気係数で非constant linear Casimirは存在しない**。
以上より両stableは、二つの相異Poisson超平面がoldhへ合体する任意の正則quadratic族の端点にならない。
これは固定frame／period profileを仮定しないwhole collision障害である。

## 4. line-conicの全factor assignments

target Dは `h²*L*Q3`、Q3はC上既約rank三、h,L相異。
source factor flagのquadric Fは必ずQ3へ移り、三marked linear factorsはh,h,Lへ移る。
sourceでは三Hが相異なるので、二つがhへ合体する。
§3が全磁気でこれを排除する。よってstable_line_conic全六源はwholeT02外。

## 5. three-linesの全factor assignments

base eigenformsは

```
L1=3h+x, L2=-3h+y, L3=-3h+x-y,
h=(L1-L2-L3)/9 = sum_i c_i Li, all c_i!=0.
```

全oldh-Casimir logfields `Xi=π#dlog(Li/h)` はfiberのみのlinear fieldである。
Cartanの二同時eigenvalueをlambda_i in C²として

```
lambda_1=(2,1), lambda_2=(-1,-2), lambda_3=(-1,1),
Xi=sum_j c_j*(lambda_i-lambda_j)*Lj · partial_q.     (2)
```

これは全magnetic Rに依らない。
相異i,j,kに対しLk=0上では

```
Xi=c_j*(lambda_i-lambda_j)*Lj · partial_q,
Xj=-c_i*(lambda_i-lambda_j)*Li · partial_q.
```

同一fiber方向でも係数Li,Ljはそのhyperplaneで独立であるから、
**非零constant係数のどの線形結合もLk上で恒等零にならない。**
base成分零の場なので、generic base点ではEulerに比例することもない。

target D=`h²*L1*L2*L3` でF0が取り得る全quadricは次だけである。

1. `F0=LiLj`。残る三marked linear factorsはh,h,Lk。§3の衝突障害で排除。
2. `F0=h²`。marked factorsは三Li。
   H=Li,Z=Lj,W=Lkとすれば、relative fieldsは
   `Y_Z=Xj-Xi, Y_F=-2Xi`。
   これらのLk上のconstant-combinationはXi,Xjの可逆な組合せなので、(CZ)が失敗する。
3. `F0=hLi`。残るmarked factorsはh,Lj,Lk。
   H=hならZ=Lj,W=Lkとして `Y_Z=Xj,Y_F=Xi` で(CZ)が失敗。
   H=Lj、Z=h,W=Lkなら `Y_Z=-Xj,Y_F=Xi-2Xj` で(CZ)が失敗。
   W=hの場合は(CW)を用い、同じ二つのfieldsの制限で失敗する。
   H=Lkの場合はj,kを交換する。

relative fieldの式は `Y_F=X_factor1+X_factor2-2X_H`、`Y_Z=X_Z-X_H` から従う。
Xiはbase formsを全て殺すので、Euler項に隠れた追加係数は無い。
quadricの平方、sourceHとの共通factor、全markedH順序を含む全assignmentを尽くした。
従ってstable_three_linesの全六magnetic sourceもwholeT02の外である。

## 6. 射程と次の親

本証明はtrue四period/edge sourceの二つのintrinsic pencil恒等式と、一般hyperplane-collision補題を使う。
固定profileの不整合を全parentへ外挿していない。
新point／operator／normal jetは無し、全zero-magnetic境界も残す。

T03ではH-Casimir sourceのcross係数がT02のopposite-pairと異なり、Y_F=aUの同じ恒等式を持たない。
したがって本T02陰性をT03へ移植しない。
T03/T16/T17/T36/T37などのactual positive whole-source入口又は未知親は未判明のままである。
