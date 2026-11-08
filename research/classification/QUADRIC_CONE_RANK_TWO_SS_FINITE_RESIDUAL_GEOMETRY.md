# Rank二quadric cone・SS：既約Dに残り得る三whole源への純還元

2026-11-08、PRIMARY E。新operator／profile／jet／CAS無し。
`QUADRIC_CONE_RANK_TWO_SS_BASE_REDUCTION.md` の全B,T正常形を出発点とする。
本稿は全既約D源の**必要な有限入口**を得る。残三源の存在又は既知帰属は未証明。

## 1. 全Ore Pfaffianとmodular恒等式

chart h=1、三base変数(a,b,c)、q=ac−b²を使う。
P=B+(Y+vT) wedge ∂v、全Yのdegree≤3で、Y3=−q*Euler3。
前noteでkappa≠−1に得た形を

```
B_ab=−ab,
B_ac=d−w*q−2b²,       w=1+2kappa,
B_bc=−bc,
T=rho*Euler3−w*H,     H=a∂a−c∂c,
rho*d=0,             rho=1+sigma+kappa
```

と書く。column divergenceはZ_B=(w+1)H。
全actual chart Pfaffianは

`Faff=G+v*K, G=(B_bc,−B_ac,B_ab) dot Y, K=(B_bc,−B_ac,B_ab) dot T`。

rank四PoissonのPfaffian/modular恒等式は `P*grad(Faff)=−Faff*div(P)`。
これはinverse volumeのLie derivative又は既受理BV恒等式から得られる。
div(P)のbase成分はZ_B−Tなので、v係数とconstant係数を比較すると

```
B#dK=−K*Z_B,
B#dG+K*Y=−G*Z_B+G*T.
```

K≠0のときg=G/Kとおけば、従って正確に `Y=g*T−B#dg`。
これはdense fieldでの恒等式を使うだけで、rational v-shearをactual PGL mapとは呼ばない。
以下ではこの式とYのpolynomial性だけを必要条件として使う。

## 2. rho≠0、kappa≠1/2の場合のpole-safe係数比較

rho≠0ならd=0で、K=rho*(2kappa−1)*b*q。
kappa≠1/2ならK≠0。
既約Fはb因子を持てないのでG|b=0≠0。
gのb=0でのsimple-pole residueをg0(a,c)とする。

`g0=G|b0/[rho*(2kappa−1)*ac]`

である。b=0上ではbvec=(0,w*ac,0)なので、これは正確に

`g0=w*Y_b|b0/[rho*(2kappa−1)]`。

従ってg0はa,cの**多項式**、degree≤3、かつ非零。
分子のtangent derivativeを捨てていない点が重要である。
実際B#db=−bHであり、Y=gT−B#dgのpole coefficient全体は

`g0*(T−H)−(B#dg0)|b0=0`。

B_ac|b0=−w*ac。g0の一つの非零monomial a^i*c^jについて二成分を比較すると

```
rho=w*(1−j)+1,
rho=w*(i−1)−1,
従って w*(i+j−2)=2.
```

i,j≥0、i+j≤3である。i+j=2は不可能。
残るsumは0,1,3で、順にw=−1,−2,2。
w=−1すなわちkappa=−1ではi=j=0からrho=0となり、この枝に反する。
w=2はkappa=1/2でK=0、別扱いである。
従ってこのwhole枝で既約Dが残り得るのは正確に

`kappa=−3/2, rho=+1又は−1, d=0`

だけである。他はGがbで割れ、Kもbで割れるのでactualFにb因子がある。
有限元の式へspecial parameterを代入して得た標本結論ではなく、全sourceの必要条件である。

## 3. rho=0：全Fのweight解

この枝ではK=0、Faff=G。
上のmodular恒等式のb成分はb*H(G)=0なのでH(G)=0。
従ってGはq,bの多項式、weighted degree≤5（qのweight二、bのweight一）。
残るa,c成分は

`[(d−w*q)∂q−b∂b]G=−(1+2w)*G`。

w≠0ならQ=q−d/wとおく。全多項式解はmonomial b^p*Q^j、

`p+w*j=1+2w, p+2j≤5, p,j≥0`

の線形結合である。
既約Dがbで割れないためにはp=0の項が必要で、j≤2よりw=−1/2又は−1だけ。

* w=−1/2、kappa=−3/4：全解はalpha+beta*b*Q²。
  固定したY3からbeta=5/2なので、global actualFは
  `alpha*h^5+(5/2)*b*(q+2d*h²)^2`。
  alpha=0はproperfactor、alpha非零を残す。
* w=−1、kappa=−1：全解はalpha*Q+beta*b*Q²で、Q因子がある。
  後述のunimodular枝のrho0にも同じ形を正当に適用出来る。

w=0では元のoperatorはd∂q−b∂b。
G=Σ b^p*g_p(q)とするとd*g_p'=(p−1)g_p。
有限degree多項式ではp≠1のg_pは零、p=1だけが残るため、やはりb因子。
d=0でも同じ結論であり、dを反転した補集合は無い。

## 4. kappa=1/2の全rhoは空

rho=0は前節で既にb因子。rho≠0ならd=0、K=0である。
全modular恒等式をmonomial a^i*b^j*c^kへ課すと

```
i−k=rho,
j+2k=5−rho,
従って i+j+k=5.
```

rhoが整数でなければG=0。
rho>0なら全項がaで、rho<0ならcで割れる。
rho=0ならi=k、j=5−2kが正の奇数なので全項がbで割れる。
従ってこの全枝もactual既約quinticでは空。

## 5. kappa=−1のunimodular全補集合

前のbase noteによりB=Jac(−bq+D*a−d*b+A*c)、T=rho*Euler3+H、
全T0零で、(1−2rho)D=rho*d=(1+2rho)A=0。
rho≠0,±1/2ならA=d=D=0なので§2が適用でき、空。
rho=0ならA=D=0で§3のw=−1が適用でき、Q因子で空。
従って残り得るのは

`rho=−1/2, A≠0, d=D=0`、又は `rho=+1/2, D≠0, d=A=0`。

指定constantが零なら§2へ戻って空。
unit constantの大きさはactual同時base scaling(a,b,c)→r(a,b,c)、v→r²vで
A又はDをr²倍へ動かして一に出来る。逆はr非零、有限平方根だけである。
全Yを同時pushし、特定unknownを零にしない。

## 6. 次に残す有限whole源

従ってconstant-C semisimpleの全既約quintic源は、全Y compatibilityを残した次だけへ入る。

1. kappa=−3/2,rho=±1、d=0の二fixed B,T。
2. kappa=−3/4,rho=0、d任意、alpha非零で上記actualF形を持つ全源。
3. kappa=−1,rho=∓1/2、対応A又はDを一とした二fixed B,T。

反転した全unitの零補集合は上でproperfactor又はG=0へ戻した。
左右reflectionの同値を未証明で使って個数を減らしていない。
次は元全sourceのY linear/Pf identities又は既知actual parentの構造認識が必要である。
本稿だけから残源の存在、componenthood、既知帰属又はE完了を主張しない。
