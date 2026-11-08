# R=3L・m二SF01：全irregular源を五つの固定P0へ還元する

2026-11-08、root純手証明、独立確認前。新CAS/operator/jet/GB無し。
actualD被約既約quintic、rank二restriction plane、f=a³、normal位数二、
generic SF normal quadratic、primitive根degrees(0,1)の全sourceを扱う。
必要first-symbol還元だけであり、五つの全higher源のG/Eは未結論。

## 1. 全closed formsと実P0係数

既pole-budget入口の全sourceはconstant normal Borelとbase座標で

```
G2=gamma*a²*u*(b*u+a*v), gamma!=0,
M=[[1,0],[b,a]],
eta=theta1-theta2=-dlog a+d(L/a), L=r*b+s*c,
beta11=theta, beta22=theta-d(L/a),
beta12=0, beta21=-d(b/a)-(b/a)*d(L/a).
```

thetaは全dlog a+d(P2/a²)、全五irregular係数を保持する。
closed diagonalとtracefree gaugeからXa=Ya=0。
`iE iXvol/f=dlog f-3beta11`と全off-diagonal contractionを直接戻すと

```
Y-X=(0,-3s*a,3r*a),
Z=(s*a,-2s*b,s*c+3*a+3*r*b).
```

実際iE iZvolの三係数は
`(-3b(a+L),3a²+3rab,3sab)`で、-3f beta21と一致する。
全normal括弧はk*u²+kappa*u*v、v²係数はroot保存とtargetから零。
全target比較前のXはa-column二係数とtracezero bc-block全三係数を持つ。

Euler-Pfaffianの全u²/uv係数を戻すと、全P0は必要に

```
X=(0,-kappa*b+p*a,C*b+kappa*c+q*a),
Y=X+(0,-3s*a,3r*a),
k=p+s,
s*kappa=0, r*kappa=s*C,
gamma=3*(kappa-r*p-s*q)!=0
```

を満たす。導出は次の全polynomial恒等式で確認できる。

```
det(E,X,Y)=3a²*(r*Xb+s*Xc),
det(E,X,Z)=3a*((a+r*b)*Xb+s*b*Xc), Z(f)=3s*f,
G2[uv]=3f*kappa-det(E,X,Y),
G2[u²]=3f*k-det(E,X,Z)-Z(f).
```

これらをgamma*a³とgamma*a²*bへ一致させると、上の全式が出る。
normal K又はtranslationを先に消したsliceではない。

## 2. kappa非零の全源は元bad-vector curlで空

s*kappa=0、r*kappa=s*Cより、kappa非零ではr=s=0、従ってY=X。
全P1三十/P2九のconstant-root bad vectorCbadは非零が必要である。
元projective curlは

`(2*kappa*I+3*DX)*Cbad=0`。

DXの三固有値は0,-kappa,kappaであるから、係数行列の固有値は
2kappa,-kappa,5kappa、全て非零。従ってCbad零となりactualDにu因子、矛盾。
C、p、q全係数を保ち、DXを最初からsemisimpleに限定していない。

以後kappa零。gamma非零からr*p+s*q非零、従って(r,s)同時零ではない。

## 3. s零の全源はnilpotent ratio invariantで空

s零ではr,pとも非零、k=p非零、C,qは任意である。
全fieldsは

```
X=(0,p*a,C*b+q*a), Y=X+(0,0,3r*a),
Z=(0,0,3*(a+r*b)), k=p.
```

全raw36/P2九とrow divergenceを戻し、canonical strict Casimir liftを用いる。
zeta=v/uとするとM=Pi sharp d(v/u)はzetaを保ち、
normal blockは3k/5倍Euler、base blockは

`A(zeta)-2kI/5+B(u,v), A=Z-zeta*(X-Y)`。

Aは`3*((1+r*zeta)*a+r*b)*partial_c`、C(zeta)上nilpotentである。
k非零なのでkI-Aは可逆、全normal degree一translation Bを消せる。
normal total degree mの全F summandへのoperatorは
Sym^(5-m)(A)+k(m-2)Iで、m二以外可逆。
全Fは元leading gamma*a²*u*(b*u+a*v)をこの有理base translationで戻したものに等しい。
C(u,v)上で二つ以上のpositive-base-degree factorsを持つため、Gauss lemmaでpolynomial Fの既約性に反する。
全C,q、全higher、全irregular係数を保った全源空性である。
有理translationをglobal PGLと呼ばない。

## 4. s非零の全源は五つのspectral壁だけ

s非零ではC=0、k=p+s、gamma=-3(r*p+s*q)非零。
同じratio fieldのA(zeta)は全zetaでtriangular、固有値(s,-2s,s)である。
base degree d上の全固有値はs*(d-3j)、j=0,...,d。
normal m=3/4/5の全operator failureは

```
m=3,d=2: k/s∈{-2,1,4},
m=4,d=1: k/s∈{-1/2,1},
m=5,d=0: k/s=0.
```

translation failure k/s=1,-2もこの和に入る。
この五値以外では前節と同じ全leading-factor/Gauss矛盾で全空。
k零をunitと仮定せず、そのwhole壁も残す。
従ってp/s∈{-3,0,3,-3/2,-1}が全必要値である。

actual base shear old c=c+t*bはr'=r+s*t、q'=q-t*pを与える。
t=-r/sでr'零、q'=q+r*p/s=(r*p+s*q)/s非零に出来る。
全normal root rowはb*u+a*vのままで、inverseはminus shear。

old base(a,b,c)=(alpha*a,beta*b,chi*c)、old normals(u,v)=(mu*u,nu*v)、
nu=beta*mu/alpha、tensor scale sigma=beta*chi/alpha²と取れば、
anchor3a²とroot rowを保ち

```
p'=chi*p/alpha, s'=chi*s/alpha,
q'=beta*q/alpha, r'=beta*r/alpha.
```

chi=alpha/s、beta=alpha/qでs'=q'=1、r'=0となる。
全unitの逆が存在し、nonzero gammaも保存する。
従って必要残余は次の五つの固定P0だけである。

```
p∈{-3,0,3,-3/2,-1},
X=(0,p*a,a), Y=(0,(p-3)*a,a),
Z=(a,-2b,c+3a), k=p+1, kappa=0,
actualG2=-3*a²*u*(b*u+a*v).
```

全higher P1三十/P2九はこのactual inverseで同時に戻す。
残五源の存在、G上界/E帰属は別の有限義務であり、今回の純必要還元から完成とは言わない。
R3L SF00 whole(t,d)、SF02既空、normaldouble既八23被覆、R=a²b、m三及び全分類は別に保持する。
