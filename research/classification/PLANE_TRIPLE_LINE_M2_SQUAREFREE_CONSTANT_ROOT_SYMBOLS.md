# R=3L・m二SF00：全irregular接続から五固定源と一つの全二parameter源へ

2026-11-08、PRIMARY G、root純手証明。独立review前、higher計算未実行。
actualD被約既約quintic、rank二restriction plane、f=a³、normal位数二、
generic squarefree normal quadraticの両primitive根がconstantの場合を扱う。
全sourceの必要first-symbol還元だけを主張し、G上界/E帰属はまだ結論しない。

## 1. 全irregular poolを落とさないraw P0

既全pole-budget入口より、normal右GL2でG2=gamma*a³*u*v、gamma非零に出来る。
全connectionはdiag(theta1,theta2)、各theta_i=dlog a+d(P_i/a²)、P_iは任意quadratic modulo C a²。
tracefree linear mixed fieldsをX,Yとして、全P0は

```
{b,c}=3a², {a,b}={a,c}=0,
{base,u}=uX, {base,v}=vY,
{u,v}=kappa*u*v.
```

theta_i closed及びtrace X=trace Y=0から、直接Cartan式
d(iE iX vol/f)=X(f)iEvol/f²によりX(f)=Y(f)=0、従ってX(a)=Y(a)=0。
これはfがa³の場合だけの結論であり、f=a²bのdiagonal fieldsへ延長しない。
normal u²/v²係数はtargetに無く、両off-diagonal mixed fieldsも零なのでnormal括弧のその二係数も零。

a-rowが全て零なのでprojective Jacobiはstrict Jacobiと同値である。
元quad(a,b,u,v)等はa倍のstrict cyclic Jacobiそのものとなる。
従ってa=1上のaffine fieldsX,Yはcommuteし、それぞれdivergence零のaffine Hamiltonian fieldsである。
全constant translationsを保持し

```
X_bc=H1*(b,c)^t+a*l1, Y_bc=H2*(b,c)^t+a*l2,
tr H1=tr H2=0,
[H1,H2]=0, H1*l2−H2*l1=0.
```

二HamiltonianのPoisson括弧はconstantなので
det(X_bc,Y_bc)=d0*a²。actualD2のuv係数から

`gamma=3*kappa−d0 != 0`。

## 2. kappa非零の全源は四つの半単純P0

全raw36/Euler六とP2九を保持する。二constant根それぞれのbad vectors C,C'は非零が必要であり、
元projective curlは

```
(2*kappa I+4DY−DX)C=0,
(-2*kappa I+4DX−DY)C'=0.
```

kappa非零なら両式のa-rowからCa=C'a=0。
commuting tracezero二行列H1,H2が共にnilpotentなら、前式は2*kappa I+nilpotentで可逆となるため不可能。
従って一方は非零半単純で、同時対角化して

`H1=diag(lambda,-lambda), H2=diag(mu,-mu)`

と出来る。zero行列の場合も以下の二条件が排除するので先にunitとは仮定しない。
bad directionsの固有値条件から、あるepsilon,epsilon'∈{1,-1}で

```
4mu−lambda=2epsilon*kappa,
4lambda−mu=2epsilon'*kappa.
```

よって(lambda,mu)は共に非零で、mu=lambdaかつkappa/lambda=±3/2、
又はmu=−lambdaかつkappa/lambda=±5/2だけである。
H1可逆とcommuting affine式から同じactual base translationでl1,l2を同時に消せる。
従ってd0零、gamma=3kappa非零。
base diagonalとtensor scaleを組み合わせ、{b,c}=3a²を保ちながらlambdaを一に出来る。
全sourceは次の四固定P0へ含まれる。

```
X=b partial_b−c partial_c,
(Y,kappa)=(X,3/2),(X,-3/2),(-X,5/2),(-X,-5/2).
```

同値な固定源がこの安全な表で重複しても被覆には差し支えない。
この段階では他の全higher方程式を解いていない。

## 3. kappa零の全源：translation固定源又は二parameter nil源

gamma非零はd0非零を要求する。
H1/H2に半単純な非零元がある場合、共通fixed pointへ移すとd0零となり矛盾する。
従って両Hはnilpotentである。

両H零ならl1,l2はindependent。base GL2とtensor scaleの実可逆変換により

`X=a partial_b, Y=a partial_c, kappa=0, gamma=−1`

という一固定P0へ戻せる。
実際translation pairを行列L、det L=delta非零とし、sigma=1/delta、base行列g=sigma Lと取る。
sigma/det g一なのでanchor3a²を保ち、sigma g^-1 L=Iとなる。

一方のH非零なら二根交換でH1非零とする。sl2のnilpotent centralizerからH2=t H1、tは任意複素数。
commuting affine式はD=l2−t*l1∈ker H1=im H1。
d0=det(l1,D)非零なのでD非零かつl1はそのline外である。
base affine shearでl1のimage方向成分を消し、base diagonal/tensor scaleのnonzero unitsを合わせると

```
X=a partial_b+b partial_c,
Y=t*(a partial_b+b partial_c)+d*a partial_c,
kappa=0, t∈C, d∈C*, gamma=−d.
```

が全sourceを覆う。具体的にXb=r*a、Xc=s*b、r,s非零へ移した後、
old b=beta*b、old c=chi*c、tensor sigma=beta*chiとし、chi=1/r、beta²=1/sを選ぶ。
両nonzero square-root choicesを許せばXの二係数は一になり、dは新しい任意非零parameterとして残る。
**dを一へ固定していない。** t任意、d非零の全二parameter源が残る。
逆変換はminus shear/nonzero diagonalの逆であり、Hの零、translation rank低下、gamma零は上記別枝に処理済み。

## 4. 次の最弱義務

全SF00は四半単純固定源、一translation固定源、一つのwhole(t,d) nil源へ必要還元された。
全P1三十/P2九は未制限である。独立幾何review後、全六源を一回の有限必要dimension/factor batchで扱うなら、
二parameter nil源の全rank/denominator壁を含める必要がある。
sample(t,d)、d一のslice、あるkernel vectorだけの因子では全源排除にならない。
今回CAS/operator/minor/jet/GBは実行していない。
他flags、R=2L+M、normal double、m三、全G/Eは本稿の結論ではない。
