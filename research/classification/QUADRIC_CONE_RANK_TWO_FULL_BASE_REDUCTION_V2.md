# Quadric-cone rank二：省略のない全Bからの正常形 V2

2026-11-08、PRIMARY E、純手証明。新operator／CAS無し。
入口は `ENTRANCE_ERRATUM` の正しい全33 sourceであり、旧非零C六係数入口を使わない。
必要supersourceの全rank四排除に使うので、actual affine変換の全逆を保持すれば十分である。
変換後の全Yを未解決のまま残し、元sourceの特定unknownを零に置かない。

## 1. SSの全JacB係数

erratumの記号A,d,D,r,s,t,xi,zeta,eta,nu、kappa、w=1+2kappaを使う。
全bvec dot curl(bvec)/2のdegree二部分は

```
s*ac+(w−2)s*b²−(kappa*t+w*zeta)*ab−(kappa*r+eta)*bc
−4(kappa+1)xi*a²−4(kappa+1)nu*c².
```

従ってs=0、(kappa+1)xi=(kappa+1)nu=0、
kappa*t+w*zeta=0、kappa*r+eta=0が全sourceで必要。

### 1.1 kappa≠−1

w非零ならxi=nu=0、eta=−kappa*r、zeta=−kappa*t/w。
actual base変換a'=a+r,c'=c−t/w,b'=bを取ると

```
B_a'b'=A−a'b',
B_a'c'=d−tr−w*a'c'+(w−2)b'²,
B_b'c'=D−b'c'.
```

w=0（kappa=−1/2）ならt=xi=nu=0、eta=r/2、zeta任意。
actual a'=a+r,c'=c−2zetaで同じcanonical形（w0,d不変）となる。
両変換は全P4のlinear map、逆は負translationで、割ったのは表示したw非零branchだけ。

column divergenceは2(kappa+1)H、H=a'∂a'−c'∂c'。
JacBの残りからA=D=0。以後primeを省きdを新constantとする。
T1=rho*Euler3−wH、rho=1+sigma+kappaは元sourceの全linear part。
Lie_T B=0とdivT constantから[T,H]=0、ab係数から全T0零。
最後にconstant coefficientがrho*d=0を与える。

従って**kappa≠−1の全sourceが**旧SS finite noteのcanonical B,Tへ戻る。
追加Fu-only係数を零としたのではなく、元JacBとactual affine逆で吸収した。
全Yのdegree≤3とleading Y3=−q*Euler3はこれらのaffine変換で保存される。

### 1.2 kappa=−1

必要式はs=0,zeta=−t,eta=rで、xi,nuを消してはいけない。
全Bは正確にJacobian(f)、

`f=−bq−t*ab−r*bc+2xi*a²−2nu*c²+D*a−d*b+A*c`。

T1=sigma*Euler3+H。Lie_TB=0はgrad(Tf−3sigma*f)=0。
degree二を全比較して

```
xi(2−sigma)=0, nu(2+sigma)=0,
T0_a=(1+sigma)r, T0_b=0, T0_c=(sigma−1)t.
```

actual a'=a+r,c'=c+tでT0を全sigmaについて消す。
fのab,bc termsも消え、xi*a²／nu*c²のtranslation由来linear項は新A,d,Dへ吸収する。
degree一は `(1−2sigma)D=sigma*d=(1+2sigma)A=0`。
したがって旧canonical cubic basesに加え、必ず次の二whole例外を保持する。

* sigma=2,xi非零：f=−bq+lambda*a²、lambda=2xi非零、T=diag(3,2,1)。
* sigma=−2,nu非零：f=−bq+lambda*c²、lambda=−2nu非零、T=diag(−1,−2,−3)。

lambda零は旧canonical cubic sourceへ戻る。lambdaを一にしなくても以下の排除に十分である。

## 2. nilの全JacB係数

erratumのxi=p26,eta=p27,zeta=p29,omega=p30を含む完全式を使う。
bvec dot curl(bvec)/2のdegree二部分は

```
[−kappa*r+(4kappa−2)zeta]*c²
+[2kappa*s+(8kappa−6)eta]*bc
+[t+(8kappa−12)xi]*b²
+[−(1+kappa)t+8kappa*xi]*ac.
```

従って

```
t=(12−8kappa)xi,
(2kappa+3)(kappa−1)xi=0,
kappa*r=(4kappa−2)zeta,
kappa*s=(3−4kappa)eta.
```

m=1−kappa、R=r−2zeta、S=s+eta、L=t−4xiとおけばcolumn divergenceは

`Z_B=2(R−2m*b,S−m*c,L)`。

### 2.1 m非零、xi=0の全source

この枝はkappa=0も含み、kappaを割らない。
actual b'=b−R/(2m),c'=c−S/mでZ_Bを−2m*N、N=2b'∂a+c'∂b'へcenterする。
全Y,Tをpushし、leading B2,Y3,T1は不変。
変換後の全linear B係数と全constantを残し、divergenceのconstant零、
JacBと同値なbvec dot N=0を比較すると、正確に

`B=B2(kappa)+ell*Jac(q)+h*c*(∂a wedge ∂b)+A*(∂a wedge ∂b)`

に限られる。ここでB2の上三角は(2kappa*ac+(1−2kappa)b²,2bc,c²)。
この表示には元Fu-onlyが残すell,hを含み、旧六係数表示と同じだと仮定していない。
実際全linear coefficientsを比較すると、B_acのlinear部は2ell*b、
B_bcはell*c、B_abはell*a+h*cだけ、constantはB_abのAだけになる。

T1=sigma*Euler3+w*N、w=1−2kappa。
[T,Z_B]=0よりT0=tau∂a。
Lie_TBの全係数は正確に

`sigma*ell=0, 2kappa*tau−sigma*h=0, ell*tau−2sigma*A=0`。

sigma非零ならell=A=0。
kappa非零ではactual a'=a+h/(2kappa)がhとtauを同時に消す。
kappa零ではh=0、actual a'=a+tau/sigmaがtauを消す。
従ってsigma非零全sourceは旧canonical B2,T1へ戻る。

sigma=0,kappa非零ならtau=0、actual a'=a+h/(2kappa)でhを消す。
新Aへconstantを吸収し、B=B2+ell*Jac(q)+A∂a wedge ∂b、T=wNを全て残す。
sigma=kappa=0ではell*tau=0、hも残す。これを次の排除V2で分岐する。

### 2.2 kappa=1（m=0）の全unimodular source

r=2zeta,s=−eta,t=4xiであり、全BはJacobian(f)、

`f=cq−4xi*ab−eta(ac+2b²)−2zeta*bc−omega*c²+D*a−d*b+A*c`。

T1=sigma*Euler3−N。
sigma非零ならT1は可逆なのでactual affine translationで全T0を消せる。
grad(Tf−3sigma*f)=0をdegreeごとに比較する。
degree二／一のoperatorsは−sigma−N、−2sigma−Nで可逆だからf2=f1=0。
従ってB=Jac(cq),T=T1へ戻る。

sigma零ならdegree二がxi=0、T0=(−2zeta,3eta,0)を強制する。
actual b'=b+zeta,c'=c−3etaでT=−N'。
Tf constantを原点で評価するとconstant零だからfはC[c',q']に入り、leading termはc'q'。
全degree≤3よりf=c'q'+ell*q'+h*c'²+A*c'+constant。
actual a-translationでhを消し、constantを捨てれば

`f=(c+ell)q+A*c, B=Jac(f), T=−N`。

この全ell,Aも消してはいけない。次の排除V2で全て扱う。

### 2.3 kappa=−3/2、xi非零の新例外

これは旧稿に無いwhole branchである。
t=24xi、L=20xi、Z_Bのc-componentは40xi非零。
T1=sigma*Euler3+4N。
[T,Z_B]=0のc成分がsigma=0を強制し、残りから

`T=−(4/5)Z_B+delta∂a`。

Z_BはPoisson modular fieldなのでLie_ZB B=0。
Lie_T B=0と∂a B_abのc係数−3からdelta=0。
従ってK=bvec dot T=0（JacBはbvec dot Z_B=0）。
全Faff=Gのmodular恒等式はB#dG=−(9/5)G Z_B。
skew性でdGを掛けると、G非零ならZ_B(G)=0。

actual b'=b−R/5,c'=c−2S/5で
`Z_B=(−10b',−5c',40xi)`。
これはhomogenization後regular nilpotent Jordan4。
全degree≤5不変polynomialの最高degree五部はc'⁵のconstant倍だけである
（次の排除V2に完全六基底と証明を記す）。
しかし元actualGの最高部は固定Y3から

`G5=2c*(ac−b²)²`。

actual affine translationでこの最高部は変わらず、a²c³係数が非零。
上の不変polynomialにはなれない。従ってこの新whole例外はactualD零又は不整合で空。

## 3. 証明のscope

全変換はactual affine linearで全P4へ延び、表示unitの逆と負translationが明示逆。
parameterを特定値へsampleせず、全Fu-only項／Y／T0を含めた必要sourceを正常化した。
残るwhole空性は別 `IRREDUCIBLE_QUINTIC_EXCLUSION_V2` のpole／ideal／valuationへ接続する。
旧稿を黙って書き換えず、erratumを保存した上での修復である。
