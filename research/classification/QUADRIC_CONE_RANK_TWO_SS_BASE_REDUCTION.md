# Rank二quadric cone・constant C semisimple：全B,Tの必要正常形

2026-11-08、PRIMARY E、純手計算による全source削減。
**以下は必要sourceの還元であり、全Y compatibility／E帰属は未解決。**
新operator、kernel、minor、profile、jet又はCASは実行していない。

## 1. 保存全33係数の完全base式

`QUADRIC_CONE_RANK_TWO_SEPARATE_REVIEW.md` の正しいOre入口を使う。
constant C≠0のsemisimple型を正規化し、a=u0,b=u1,c=u2、q=ac−b²とする。
fullsourceは `Pi=B+(Y+vT) wedge partial_v`、degY≤3、degT≤1。
保存PENCIL report cases[0]の全Bを、前のC零noteと同じ六係数

`A=p1,d=p2,D=p8,r=p7,s=p14,t=p20`

及びFvの係数kappa=p31で読むと

```
B_ab=A−r*b+s*a−a*b,
B_ac=d−r*c+t*a−(1+2kappa)*a*c+(2kappa−1)*b²,
B_bc=D−s*c+t*b−b*c.
```

他の全family係数はY,Tへ残し、零にしない。
Tの全linear partは保存sourceから

`T1=diag(sigma−kappa,1+sigma+kappa,2+sigma+3kappa)`。

T0はまだ全三constant。以後の操作は元fullJacobiの必要条件JacB=0、Lie_T B=0だけを使う。

## 2. JacBの全係数による三branch

bvec=(B_bc,−B_ac,B_ab)として、bvec dot curl(bvec)/2は

```
s*a*c+(2kappa−1)*s*b²−kappa*t*a*b−kappa*r*b*c
−(1+kappa)*D*a+(1+kappa)*A*c+s*d−r*D−t*A.
```

従って常にs=0。

* kappa≠0,−1：r=t=A=D=0。
* kappa=0：A=D=0、r,tを全て残す。
* kappa=−1：r=t=0、A,d,Dを全て残す。

この全splitにparameter unitを隠していない。
kappa=0でactual affine base変換a'=a+r,c'=c−t,b'=bを行うと

`B_ab=−a'b, B_ac=(d−rt)−a'c'−b², B_bc=−bc'`。

これは全P4のlinear PGL変換(a'=a+rh,c'=c−th)で、逆は同じtranslationの負。
Y,Tは全体をpushし、未決のhigher係数を固定しない。
Sはh=0で不変、leading qとleading T1は変わらない。
変換後の下位Y係数を既定sourceの特定成分へ無断で同定したものではない。

## 3. kappa≠−1の全必要B,T

前節のactual affine変換後、kappa≠−1の全Bは

```
B_ab=−ab,
B_ac=d−(1+2kappa)ac+(2kappa−1)b²,
B_bc=−bc.
```

column divergenceは `2(1+kappa)*(a∂a−c∂c)`。
TがaffineでdivTはconstantなので、Lie_T B=0とBV恒等式からTはこのlinear fieldと可換。
そのweight−1,0,1は相異なるため、T0のa,c成分は零。
T0のb成分をdeltaとすると、Lie_T Bのab成分は−delta*aだからdeltaも零。
従ってT=T1である。

rho=1+sigma+kappa、H=a∂a−c∂cと書くと

`T=rho*Euler3−(1+2kappa)*H`。

Lie_T Bのac成分のconstant coefficientは−2rho*d。
他のhomogeneous coefficientは恒等的に零なので、必要十分条件は

`rho*d=0`。

kappa=1/2のb²係数零も含む全kappa≠−1で成立し、そこで割っていない。
したがって残る全sourceはd=0又はrho=0の二closed源に入る。
この段階では全Yを保持し、actualDのproper factorや既知parentをまだ主張しない。

## 4. kappa=−1のunimodular cubic base

全r=t=s=0のもとで

```
f=−abc+b³+D*a−d*b+A*c,
B=Jacobian(f),
T1=sigma*Euler3+H.
```

ここでJacobian(f)は(B_bc,−B_ac,B_ab)=gradfという規約。
Lie_T B=0は `grad(Tf−(divT)f)=0` と同値。
T0の三成分はfのhomogeneous cubicへの作用としてbc,ac−3b²,abの三独立quadricsを作るので全零。
残るlinear coefficientsから正確に

`(1−2sigma)*D=0, sigma*d=0, (1+2sigma)*A=0`

を得る。sigma=0,±1/2を捨てず、全三constantを保持したunit／zero分岐である。
sigmaがこれら以外ならA=d=D=0でf=−bq。

## 5. 次の最弱義務

本稿はsemisimple constant-Cの全base Jacobiと全Tを、小さいwhole正常形へ縮めた。
Y3=−q*Euler3と全Y2,Y1,Y0はまだ残り、
`Lie_Y B=Y wedge T`の全元compatibilityを満たす必要がある。
generic stage kernelをshearで除けるという推測だけで全特殊源を消していない。
nilpotent constant-C sourceは別の未処理branch。
固定のB,T正常形に対する全Yの因子／actual既知parent／有限orbit認識を次の構造的義務とする。
