# Rank二quadric coneの全既約quintic源は空

2026-11-08、PRIMARY E、producer純証明。独立reviewは別receiptに固定する。
新CAS、M1、operator、kernel、minor、profile又はjetは一度も実行していない。

**定理。** P4上のrank四Poisson構造の実退化因子Dが被約既約quinticなら、
その非正規曲面支持に、一般restriction rankが二のintegral quadric coneは存在しない。
全R multiplicity、全vertex、全higher、全normalizer zero/unit補集合を含む。
可約actualDの所属を分類したという意味ではない。

## 1. 全sourceの有限必要入口

全元150係数からのcomplete33源と、正しいHamiltonian pencil D=0の帰結は
既受理 `QUADRIC_CONE_RANK_TWO_SEPARATE_REVIEW.md` に従う。
新pure chainは次の四note。

* `ZERO_C_ACTUAL_FACTOR`：constant C零は全Y/T込みでactual quadratic factor。
* `SS_BASE_REDUCTION` と `SS_FINITE_RESIDUAL_GEOMETRY`：C semisimple非零の全源は下の三typeへ入る。
* `NIL_BASE_REDUCTION` と `NIL_FINITE_RESIDUAL_GEOMETRY`：C nilpotent非零は下の唯一のtypeへ入る。

ここでa=u0,b=u1,c=u2、q=ac−b²、h=1、P=B+(Y+vT) wedge ∂v、全degY≤3。
Yを選んだsubfamilyへ制限していない。全actualFaffはG+vKで
G=bvec dot Y、K=bvec dot T、bvec=(B_bc,−B_ac,B_ab)。
K非零では全modular/Pf恒等式から `Y=gT−B#dg, g=G/K`。

残るsourceは

1. SS kappa=−3/2、rho=±1、d=0。
2. SS kappa=−3/4、rho=0、d任意、G=alpha+(5/2)b(q+2d)²、alpha非零。
3. SS kappa=−1、rho=−1/2,A=1、又はrho=+1/2,D=1の二unimodular cubic bases。
4. nil kappa=3/4、sigma=0、A任意、G=alpha−(5/2)c(q+2A)²、alpha非零。

以下でこの全てが空となる。左右同値を証明せず一方を捨てることはしない。

## 2. SS／nilのalpha族はbase zero idealだけで空

**SS type2。** bvecは

`(−bc,−d−ac/2+5b²/2,−ab)`。

proper ideal I=(b,q+2d)で全成分が零となる。
実際b=0、ac=−2dを入れれば中央成分も零。
Iは全dでproper：商はC[a,c]/(ac+2d)で非零である。
G=bvec dot YはIに属する一方、表示Gをmod Iへ送るとalphaとなる。
従ってalpha=0、必要unitに矛盾する。

**nil type4。** bvecは

`(c²,−2bc,A+3ac/2−b²/2)`。

全AでproperなI=(c,b²−2A)に全成分が属する。
商はC[a,b]/(b²−2A)で、A=0でも非零。
G mod Iはやはりalphaなのでalpha=0、矛盾する。
有限parameterを選んだpoint検査でなく、全parameter sourceのideal必要条件である。

## 3. SS type1はq-poleの負Euler weightで空

w=1+2kappa=−2、rho=±1。
全B,Tは

```
B_ab=−ab, B_ac=2ac−4b², B_bc=−bc,
T=rho*Euler3+2H, H=a∂a−c∂c,
K=−4rho*b*q.
```

Fにq因子がないと仮定する。g=G/Kのq=0におけるsimple-pole coefficientをetaとする。
q=0ではbvec=−b*gradqなので

`eta=G|q0/(−4rho*b)=(gradq dot Y)|q0/(4rho)`。

従ってetaはgraded integral ring C[a,b,c]/(q)のpolynomial class、非零でdegree≤4。
bを反転して定義した有理式に留まらず、分子bが正確に消えることを使っている。

B#dq=2qH、B|q0=Euler3 wedge Hなので、Y=gT−B#dgのq-pole全体を消す必要条件は

`eta*(T+2H)−(B#d eta)|q0=0`。

Euler3,Hはquadricのdense smooth locus上で独立。
B#deta=Euler3*H(eta)−H*Euler3(eta)より

`H(eta)=rho*eta, Euler3(eta)=−4eta`。

graded polynomial ringのEuler eigenvaluesは非負整数で、負四のclassは零。
eta非零に矛盾する。従ってGはqで割れ、Kもqで割れるので全actualFにproper q factor。
全Y、rhoの両符号、quadric vertexもpolynomial恒等式として含まれる。

## 4. 二unimodular cubic basesは非平方unitの対数微分で空

**A-unit case。** rho=−1/2、

`f=−bq+c=c(1−ab)+b³, B=Jac(f), T=rho*Euler3+H`。

全Pfaffianのv係数は正確にK=T(f)=3rho*f。
fはprimitive c-linearでirreducibleなので、その商のfunction fieldはC(a,b)。
Fにf因子がなければeta=G/(3rho) mod fは非零である。
B#df=0なのでYのf-poleを消す必要条件は

`B#deta=eta*T` on f=0。

このsurfaceのfunction fieldの座標a,bで

`B_ab=1−ab, T_a=a/2, T_b=−b/2`。

従って

```
(1−ab)*partial_b eta=(a/2)*eta,
(1−ab)*partial_a eta=(b/2)*eta.
```

これらから `partial_a[eta²(1−ab)]=partial_b[eta²(1−ab)]=0`、
ゆえにeta²(1−ab)は非零constant。
しかしC(a,b)でirreducible prime1−abのvaluationは一であり、squareのvaluationは偶数。
この恒等式は不可能で、eta=0でなければならない。
従ってGもKもfで割れる。

1−abがaffine surface上でunitであることはこの証明の障害ではない。
valuationはそのfunction field C(a,b)で取り、境界divisorも含む。
unitを無断でsquare rootしたり、surface内の零divisorを仮定したりしていない。

**D-unit case。** rho=+1/2、f=−bq+a=a(1−bc)+b³。
field C(b,c)、B_bc=1−bc、T_b=b/2、T_c=−c/2で同じ二式が成立する。
eta²(1−bc)=nonzero constantは同じodd valuationで不可能。
このcaseも元のwhole sourceのままproper cubic factorとなる。

## 5. Global actual因子と全pointwise Eへのscope

得たq factor、f factorはchart h=1でそれぞれdegree二／三。
actualFの元degreeは五なので、polynomial cofactorのdegreeは高々三／二であり、
homogenizationすると元のEuler-Pfaffianにもproper factorがある。
cofactor零はactualD零なので対象rank四でない。
Euler liftの変更はactualE wedge Pi²を変えない。

全C零／semisimple／nilpotent正規化、反転した係数の零補集合、全kappa/sigma壁を
上のcomplete whole-source chainへ戻したため、指定するrank二quadric-cone枝は全空となる。
rank零quadric-coneは別の既受理whole-factor theoremに従う。
一般βのG排除を特殊点へ外挿したのでも、非閉鎖次元testをE被覆と呼んだのでもない。
新component数は零。独立全chain読解が終わるまでglobal登録はrootのみが判断する。
