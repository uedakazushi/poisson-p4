# Rank二quadric cone・constant C nilpotent：全base Jacobiの必要還元

2026-11-08、PRIMARY E、純手計算。全Yを未決のまま保持し、E帰属はまだ主張しない。
新計算／operator／minor／profile／jetは無い。

## 1. 全元源の式

正しいOre入口と元の全33 sourceは既独立 `QUADRIC_CONE_RANK_TWO_SEPARATE_REVIEW.md` を使う。
`PENCIL` report cases[1]（type=nilpotent）の全Bを読み、a=u0,b=u1,c=u2、q=ac−b²、
A=p1,d=p2,D=p8,r=p7,s=p14,t=p20、kappa=p33とおく。
highest-symbolの必要条件でFv=kappa*cだからp28=p31=0であり、全Bは正確に

```
B_ab=A−r*b+s*a+2kappa*a*c+(1−2kappa)*b²,
B_ac=d−r*c+t*a+2b*c,
B_bc=D−s*c+t*b+c².
```

全Tのlinear partは `T1=sigma*Euler3+(1−2kappa)*N`、
`N=2b∂a+c∂b`、T0は全三constant。
全Y3=−q*Euler3と全Y2,Y1,Y0を残す。

## 2. 全JacBと全kappa零補集合

bvec=(B_bc,−B_ac,B_ab)とすると `bvec dot curl(bvec)/2` は

```
−kappa*r*c²+2kappa*s*b*c+t*b²−(1+kappa)*t*a*c
+2(1−kappa)*D*b−(1−kappa)*d*c−r*D+s*d−t*A.
```

従って常にt=0。残りは全て次の三つに尽くされる。

* kappa≠0,1：r=s=d=D=0、Aを全て残す。
* kappa=0：d=D=0、r,s,Aを全て残す。
* kappa=1：r=s=0、A,d,Dを全て残す。

この分岐はhighest coefficientの必要条件であり、kappaをgeneric値に制限したものではない。

## 3. kappa≠0,1の全T

この枝ではBは上のhomogeneous quadratic部分とconstant A*∂a wedge ∂bだけ。
column divergenceは `−2(1−kappa)*N`。
Lie_T B=0とdivT constantから[T,N]=0。
T1は既にNと可換なので、T0はNのkernel、すなわちtau*∂aに限る。
Lie_{tau∂a}Bのab係数に2kappa*tau*cが現れ、kappa非零なのでtau=0。
残るconstant係数は−2sigma*A。
従って全sourceの必要条件は

`T=T1, sigma*A=0`。

まだYを除去せず、このbase形だけからactualDの因子を推論しない。

## 4. kappa=1の全unimodular cubic base

この枝は

`f=c*q+D*a−d*b+A*c, B=Jacobian(f), T1=sigma*Euler3−N`。

Lie_T B=0はgrad(Tf−divT*f)=0と同値。
T0のf3=cqへの作用はc²、−2bc、2ac−b²の独立三quadricsを作るからT0=0。
残るlinear coefficientsを全て比較すると

`sigma*D=0, sigma*d−D=0, d−2sigma*A=0`。

sigma≠0ではA=d=D=0、sigma=0ではd=D=0でA任意。
従ってこの全枝ではf=c(q+A)、ただしsigma*A=0。
fのproper factorというだけでは四変数actualDのfactorを結論出来ないため、
全Y compatibilityは別の未完了義務として残す。

## 5. kappa=0に残すsource

t=d=D=0だけを入れた全(A,r,s)と全T0,Yを保持する。
JacBだけではr,sは消えない。
affine modular fieldをcenterへ移してもBにはunimodular linear部分が残り得るので、
nilpotent fieldのzeroを直しただけでwhole tensorの正常化が終わったとはしない。
次はこのwhole sourceの実normalizer逆、Lie_T B及び全Y条件を扱う。
特定point、非零A、semisimple closure又はconstant-rank bundleを未証明で仮定しない。
