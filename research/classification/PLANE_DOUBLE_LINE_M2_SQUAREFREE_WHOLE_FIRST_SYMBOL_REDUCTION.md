# R=2L+M・m二SF：01/02全空と00一模型への全源還元

2026-11-08、root純手証明、独立確認前。新CAS/operator/jet/GB無し。
actualD被約既約quintic、rank二restriction plane、f=a²b、normal位数二、generic SFを仮定する。
既03/12/11全空と全pole-budget入口を使い、残flags00/01/02を扱う。
結論は01/02全空、00の全first-symbolが一つの固定P0へ必要還元されることである。
00のfullhigher G/Eは未結論。

## 1. constant normal root一つの全bad-vector必要条件

H=a∂a−2b∂b+c∂cとし、全closed diagonal charactersをtracefree fieldsへ戻す。
X(f)=Y(f)=0から、全fieldsは

```
X=lambda*H+(p*a+q*b)∂c,
Y=mu*H+(r*a+s*b)∂c.
```

Hの項は全irregular d(c/a)を含み、初めからnil fieldsへ限らない。
rootuを保つ全P0はmixed uXとuZ+vY、normal K=k*u²+kappa*u*v。
uv targetはgamma*f、gamma非零なのでactualPfの全係数から

```
lambda*s=mu*q,
gamma=3*(kappa+lambda*r-mu*p)!=0.             (1)
```

を得る。Z,kはこの係数に寄与しない。

全raw36/Euler六、全P2九でbad u-row v²Cは非零が必要。
元projective全式のvector恒等式は

`E∧((DY+kappa I)C)−(2Y−X)∧C+P_(C(f))=0`。

curlと全ab/ac/bc係数を比べる。例えばab成分は

`a*(kappa+lambda-4mu)*Cb-b*(kappa-2lambda+5mu)*Ca=0`。

Ca非零ならacのc係数と合わせlambda=2mu,kappa=−mu。
Ca=Cb=0,Cc非零でもac/bcの二係数から同じ結論。
Ca=0,Cb非零ならab及びbcのa係数からkappa=4mu−lambda、p=2r+2。
従って全sourceは次の必要unionへ入る。

```
I : lambda=2mu, kappa=−mu,
II: kappa=4mu−lambda, p=2r+2.
```

unit gammaはIではmu*(2r−p−1)非零、IIでは(lambda−2mu)*(r−1)非零を要求する。
ゼロCa/Cb、mu零、lambda零を先に捨てた導出ではない。

## 2. flag00：両根の必要unionから唯一の全P0

00ではZ=k=0。v-rootにも同じ全vector式を使い、
X/Y交換とkappa→−kappaを行う。
I/I、I/II、II/Iの交差は上のunit gammaと矛盾する。
II/IIだけが残り

```
mu=−lambda, kappa=−5lambda,
p=r=−2, q=−s, gamma=−27lambda!=0.
```

actual c-shearでq,sを同時に消せる。
old c=c+t*bならq'=q+3lambda*t、s'=s−3lambda*t。
t=−q/(3lambda)を選び、逆はminus shear。
次にold c=c/lambda、tensor scale1/lambda、a,b固定でanchor P_fを保ちlambdaを一にする。
p,rはa∂cとtensor scaleの相殺により−2のままである。

従って全00源は

```
Pbase=P_(a²b),
X=H−2a∂c, Y=−H−2a∂c, Z=0,
K=−5u*v, actualG2=−27a²buv
```

という一fixed P0のfullhigher fiberへ入る。
両unit normalizer及び全inverseが存在し、higher coefficientsも同時に戻す。
一fixed pointの所属からfullfiberを推測したのではない。

## 3. flag01の全sourceは空

全M=[[1,0],[A,B]]、G2=gamma*(f/B)*u*(Au+Bv)、A一次primitive。
B=a又はb、actual Aut(R)とnormal BorelでAは残component又はcへ戻せる。
全closed formsの入口を使い、off-fieldは
`E cross Z=−3(f/B)*(dA+A*eta)`、eta=theta1−theta2。
ここでcrossの三係数は `(bZc−cZb,cZa−aZc,aZb−bZa)` と明示定義する。
収縮の合成順による符号の曖昧さを避け、以下の全literal Zはこのcross式を使う。

### B=a

eta=rho*dlog a−(1+rho)*dlog b。
diagonal差からlambda=mu=ell、p−r=3(1+rho)、q=s。
Iはunit gammaに反する。IIと(1)は

`kappa=3ell, p=6rho+4, r=3rho+1, gamma=−9ell*rho!=0`。

A=bではZ=−3rho*b∂c、actualu²比較はk=0。
zeta=v/u上のcanonical ratio fieldはbase A(zeta)=Z−zeta(X−Y)がnilpotent、
normal degree係数k_eff=k+kappa*zeta=3ell*zetaはC(zeta)のunit。
全higherを保ったnormal-degree-one translationが存在し、m≠2の全operator
Sym^(5−m)(A)+k_eff(m−2)Iは可逆。
従って全actualFは元leadingの有理base translationに等しく、C(u,v)上でpositive-degree factorsを持つ。
Gauss lemmaからpolynomial Fは既約になれず、全source空。
このratio論法はkappaを落とさずk_effを用い、全P2九/全higherを含む。

A=cではZ=diag(rho+2,rho−1,−2rho−1)。
actualu²のabc係数は−9ell(rho+1)、target係数はgamma=−9ell*rho。
ell非零なので一致せず、全P0空。

### B=b

eta=−dlog a+d((rho*b+sigma*c)/a)。全diagonal差は

```
lambda=mu−sigma, p=r−3, q=s−3rho,
sigma*s=3mu*rho.
```

最後のsはYのb∂c係数である。
Iではsigma=−mu、s=−3rho、q=−6rho、mu非零、
IIではr=−5,p=−8,kappa=3mu+sigma、gamma=18(mu+sigma)非零。

A=aではZ=(0,−3sigma*a,3rho*a)。全actualu²比較は

```
k+3lambda*rho−q*sigma=0,
gamma=3sigma*(1−p).
```

Iではk=0、右辺gamma=3mu(r−4)、(1)ではgamma=3mu(r+2)なので矛盾。
IIではflat relationから3lambda*rho−q*sigma零、ゆえk=0。
gamma比較はsigma=2mu、従ってmu非零、kappa=5mu。
ratioのbase A(zeta)は全zetaでtriangular、固有値sigma*zeta,-2sigma*zeta,sigma*zeta。
k_eff=kappa*zetaなので、その比k_eff/(sigma*zeta)=5/2。
translation又はm=3/4/5 operatorの全失敗比は{-2,1,4,-1/2,0}だけで、5/2は含まれない。
従って元leadingがpositive-degree factorsを保つ有理translationへ戻り、全higher既約性に矛盾する。

A=cではZ=(−rho*a,−3a−rho*b−3sigma*c,2rho*c)。
flat relationより3lambda*rho−q*sigma零で、全u²係数は

`3a²b(k+rho−q)+3(1−p)a³+3(1−p)sigma*a²c`。

target gamma*a²cとのa³比較はp=1、続いてa²c比較はgamma=0となり矛盾する。
I/IIの区別以前にこの全A=c源は空。
これでB=a,b双方、全primitive Aを尽くした。

## 4. flag02も全P0空

全M=[[1,0],[A2,B2]]でA2はprimitive quadratic。
全pole-budget表のB2はa²又はab、etaはそれぞれ−2dlog b又は−2dlog a。
両場合lambda=mu=ell、q=s。
Iはunit gammaに反するのでIIを使う。
B2=a²ではp=10,r=4,kappa=3ell,gamma=−9ell、B2=abではp=−8,r=−5,kappa=3ell,gamma=18ell。
いずれもell非零。actual c-shearでq=s零、base c-scale/tensor scaleでell一に出来る。
全A2はこの変換後も任意primitive quadraticとして残る。

B2=a²ではh=b、Z=3(A2_c∂a−A2_a∂c)、
beta11=4dlog a−3dlog b−d(c/a)。
恒等式G2[u²]=3f*(k−beta11(Z))とtarget−9b*A2の比較は

`(a∂a+(4a+c)∂c−1)A2=(k/3)*a²`。

全六quadratic monomialsを比べるとA2∈span(a²,ab)だけ。
従ってaがA2を割り、gcd(A2,a²)=1に反する。

B2=abではh=a、Z=3(−A2_c∂b+A2_b∂c)、
beta11=−2dlog a+3dlog b−d(c/a)。
target18a*A2との比較は

`(3a∂c+b∂b−2)A2=−(k/3)*ab`。

全六monomialsからA2∈span(ab,b²)だけ。
bがA2を割るのでgcd(A2,ab)=1に反する。
非primitive境界を捨てず、元UFDで別低degree flagへ戻す範囲と区別している。

## 5. 全scope

全R2L/m二SF01/02は空、00は上の唯一fullfiberへ必要還元された。
既03/12/11排除と合わせれば、nonreduced R2L/m二SFで残るG/E義務はその一fixed P0の全higherだけ。
現在はproducer純証明で、独立全case/sign/unit監査を受ける前に全G/E完了と呼ばない。
R2L constantdouble、m三及び全Poisson分類は別義務である。
