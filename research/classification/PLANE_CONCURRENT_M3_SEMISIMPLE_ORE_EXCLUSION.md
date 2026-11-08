# K0全SS source：actualDのweightとglobal quadratic liftによる空性

2026-11-08、PRIMARY G、producer純証明・独立review前。
generic_exhaustionと全module chartを照合し、rootも最終weighted exceptionを独立導出した。
正式な別作者receiptは別に区別する。
新CAS、first-source probe、matrix、GB又はjetは無い。
fのrootsがa,bであるという未証明入口を一切使わず、fullactualDから全SS源を排除する。

## 1. 全SSと全magneticのactual factorization

全sourceはt非零、s任意について

```
U=(u²,tuv), V=(uv,tv²),
A1=(t+s)v, B1=−(1+s)u,
Qab=aA1+bB1+C2, Qan=U,Qbn=V,Quv=0,
P=Q+∂c∧W(y), E_y∧L_WQ=0.
```

全quadraticW、元normalweight/normalorder三の条件を保つ。
C2=c0u²+c1uv+c2v²とし

```
D0=u∂a+v∂b,
S=−(t+s)a∂a−(1+s)b∂b+u∂u+tv∂v
  −c2v∂a+(c0u+c1v)∂b.
```

元fullQは正確にD0∧S。
先行whole central正常形は、三lines

`s+2=0, s+t+1=0, s+2t=0`

の外でC2零、各lineからcommon center(t1,s−2)を除けば対応する一monomialの0/1 residue。
common centerでは全C2三coefficientsを保持する。
これらのtrue shift/scaling inverseを元sourceへ戻す。

## 2. FullmoduleからactualDの二つの必要恒等式

chartu1、x=a,h=b−ra,r=vとする。
これはchartのpolynomial automorphism、Jacobian一である。
D0=∂x、projectiveQ及びc-linear vectorは

```
pi=∂x∧V,
V=(μh+C2(1,r))∂h+νr∂r,
μ=−s−2, ν=t−1,
Z=∂x,
P_chart=pi+∂c∧(Y+c∂x).
```

元全Jacobiは `L_Ypi+Y∧∂x=0`。
J=det(V,Y_normal)、Y_normal=(Y_h,Y_r)と置く。
∂xを含まないmodule係数からJ_x0。
残normal equationは

`[Y_normal,V]=(∂xY_x)V+Y_normal`。

detをVで微分すると

`V(J)=(divV−1)J=(μ+ν−1)J`。                         (1)

actualDはc非依存（PfQ零）で、chart値はsignを除きJ。
したがって(1)はfullactualD自身の恒等式。
全higherYを選んだpoint又はchosen firstjetへ固定していない。

## 3. Nonzero centralの全三walls

t非一では各magnetic wallはμ=jν、j=0,1,2で

`V=μh∂h+νr∂r+γr^j∂h`, γ0/1。

γ1ならdiagonal partV0とnilpotentN=r^j∂hは可換。
全有限polynomialJのV0-weight分解で、weightμ+ν−1以外は
nonzero scalar plus nilpotentが可逆なので消える。
残weightではNJ0、従ってJ_h0。
J_x0も戻すとJはrだけのpolynomial。
quintic homogenization D=u⁵J(v/u)はu,vのbinary quinticなので複素体上でreducible。
zeroJもrank四対象外。
全nonzero magnetic complementsを閉じた。

## 4. Scalar t一とcommon center

t1,s≠−2ではC2零、ν0、μ非零。
(1)はμhJ_h=(μ−1)J。
normalorder三でbasedegree≤二なのでh-degreejは0,1,2。
μ(j−1)=−1によりj一は不可能。
j零ならDはbinary(u,v)、j二ならDは(bu−av)² times binary linear。
どちらもproper factor又はzeroである。

common center(t1,s−2)では元Qを直接戻すと

`Q=−E_y∧(u∂a+v∂b)+C2∂a∧∂b`。

actualDはC2 times cubic、C2零ならD零。
全三magnetic parametersを保持して空である。

## 5. C2零、t非一：全weighted monomialsを尽くす

ν非零、V=μh∂h+νr∂rとする。
Jはpolynomial in h,rで、normalorder三からh-degreej≤二。
actual homogeneous quinticへ戻すと

`D=Σ c_jk*(bu−av)^j*v^k*u^(5−2j−k)`,
`j=0,1,2, k≥0, 2j+k≤5`。                            (2)

uのnegative powersを隠していない。
もし最小u-exponentが負なら、u0でbu−av=−avを代入したpole coefficientは
異なるa^j v^(j+k) monomialsの和となり、全て零でなければregularにならない。
従って(2)がglobalDの全monomialsを尽くす。

normalorderが正確に三なのでj二のnonzero termが必要。
k零ならs=−t、k一ならs=−1である。
t非一では両方が同時には起きない。

### s=−t

q=(t−2)/(t−1)と書くと、possible additional monomialsは

* j二,k零：u w²、w=bu−av（必須）。
* j一,k=q：u^(3−q)v^q w、q整数0..3のときだけ。
* j零,k=2q：u^(5−2q)v^(2q)、2q整数0..5のときだけ。

Dがwで割れないためにはj零が必要。
Dがuで割れないためには、j零がu-exponent零、又はj一がu-exponent零でなければならない。
後者はq3だが、そのときj零k6は許されずw因子となる。
従って唯一のpossible irreducible caseは2q5、すなわち

`(t,s)=(1/3,−1/3)`, `D=A u w²+B v⁵`, AB非零。        (3)

j一はq5/2で存在しない。

### s=−1

同じ分類をu/v、a/bのactual swapで戻すと唯一のpossible caseは

`(t,s)=(3,−1)`, `D=A v w²+B u⁵`, AB非零。

swap後の全tensorを1/t倍すると、SS parametersは正確に(t,s)→(1/t,s/t)。
両caseはこのactual swapとnonzero tensor scaleで同じ判定に尽くされる。
全他parameterはu/v/w proper factor又はzeroで閉じた。
generic residueからfのroot位置を推測したものではない。

## 6. 最後のbinomialはglobal quadratic Wへlift出来ない

(3)でu1を取り

`μ=−5/3, ν=−2/3, J=A h²+B r⁵`。

x=h0でJ=−ν r Y_hだから

`Y_h(0,0,r)=(3/2)B r⁴`。

ところが全global quadratic Wからのchart formulaは

```
Y_a=W_a−aW_u,
Y_b=W_b−bW_u,
Y_r=W_v−rW_u,
Y_h=Y_b−rY_a−aY_r.
```

x=a0,h0ならb0なので、正確に

`Y_h(0,0,r)=W_b(0,0,1,r)−r W_a(0,0,1,r)`。

右辺はr-degree高々三。
B非零なら上のdegree四に一致せず、B零ならDがu w²でreducible。
従ってbinomialcase全体はactual fullsourceを持たない。
swapcaseも同じ。

これはnormalization-flatness又はchosen firstsource rankではなく、
元の全quadratic tensorのglobal polynomial degree obstructionである。

## 7. 結論

全t非零、全s、全central rankwalls/common center、全quadraticWについて、
reduced irreducible actual quintic・normalorder三の**全SS源は空**。
fのroota/b仮定や未確認のfirstjet monomial classificationを依存にしない。
Jordan whole sourceはgeneric_exhaustionの別pureproofへ分担し、別に独立受理を区別する。
本稿だけで全concurrent G/E又はcomponent数を更新しない。
