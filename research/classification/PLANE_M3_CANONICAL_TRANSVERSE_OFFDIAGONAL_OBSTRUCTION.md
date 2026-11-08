# m三constant-root：canonical strict liftからの横断offdiagonal obstruction

2026-11-08。producer純証明、plane_collision_reviewによる別作者の全係数手確認を受領。
全projective Poisson tensorのcanonical liftを用い、局所のstabilizerがglobalへ延長するとは仮定しない。
以下は被約R/logarithmic性を必要としない。

## 1. 全sourceと仮定

homogeneous coordinates(a,b,c,u,v)でrank二planeをu=v零とする。
全quadratic representativeをnormalweight零/一/二に分ける。
weight零はbase gradient lift P_fと

```
P0^{i,u}=uX_i,
P0^{i,v}=uZ_i+vY_i,
P0^{u,v}=k u²
```

で、X,Y,Zはtracefree linear fields。
κuv項は零とする。X,Yの像は同一定数方向C*∂cに含まれる。
constant normal root uの全weight一を保ち、そのu-rowのu非可除項をv²Cと書く。
仮定はC=λ∂c、λ非零である。
fullP1の他全constant normal vectors及びbase bivectors、fullP2を固定しない。

結論は

`(DZ C)_a=(DZ C)_b=0`。

従ってZのa/b成分がcへ非自明に依存すれば、その全sourceは既約quinticを持てない。
最後の空性にはC非零が既約Dの必要条件であるというconstant-root因子補題を用いる。

## 2. canonical liftは全actual tensorで正当に使える

δPのj成分をΣi∂iP^{ij}とする。
quadratic bivector PではδPはtrace零linear field、δ²P零。
Eは全五変数Euler field。

`Pbar=P−(1/5)E∧δP`

は同じprojective tensorを表し、δPbar零。
E∧Vのdivergenceはlinear Vに対し5V−E divVなのでこの式は直接確認出来る。
Πがprojective PoissonならE∧[Pbar,Pbar]零。
BV identityからδ[Pbar,Pbar]零。
B=[Pbar,Pbar]についてδ(E∧B)=5B−E∧δBなのでB零、すなわちPbarはstrict Poissonである。
この論理は任意の非canonical raw liftをstrictとみなすものではない。

δとE∧はnormalweightを保つため、全Pbarのweight一Jacobiもstrictに零である。
weightを上げたformal approximationではなく、元quadratic tensorの正確な係数を使う。

## 3. weight零とbad Cはどう変わるか

tracefree X,Y,Z、κ零から

`δP0_base=−(X+Y), δP0_u=0, δP0_v=2k u`。

よってcanonical symbolは

```
Xbar=(4X−Y)/5, Ybar=(4Y−X)/5,
Zbar=Z−(2k/5)Ebase, Kbar=(3k/5)u².
```

Xbar/Ybarの像もc方向である。
Zbarがtracefreeとは仮定しない。base bivectorがgradientのままとも仮定しない。
以下で使うstrict componentはbase-base bivectorへ依存しない。

weight一のmixed normal-column係数はbaseに定数なので、δP1のnormal componentsは零。
そのEuler correctionはu-rowへu倍の項だけを加える。
従ってbad vectorは正確にCbar=Cである。
他のP1 constant vectorsは変化しても、引き続き全constant vectorsとして保持する。

## 4. strict uv² coefficient

一般の全mixed columnsを

```
g=uX+u²A+uvB+v²C,
h=uZ+vY+u²D+uvF+v²G
```

とし、K=k u²+κuvとする。
strict Jacobi component J^{i,u,v}のuv²係数は直接の微分で

`(2k I+DZ)C+(DY)B−(κI+DX)G`。

base P1/P2やA,D,Fはこのnormal monomialへ寄与しない。
canonical tensorへこの式を適用する。
DYbar B及びDXbar Gはc方向、kbar Cもc方向、κbar零。
従ってa/b成分は(DZbar C)のa/b成分だけであり、strict Jacobi零からそれは零。
Zbar−Zはscalar Eulerなので、そのC微分もc方向。
これで主張が従う。

## 5. R2L generic tripleの全r八dimへの適用

x=b/a,y=c/aとし、全rを

```
r=(Pdx+Qdy)/x,
P=A+Bx+Cy+Dxy+Ey²,
Q=F+Gx+Hy−Dx²−Exy
```

で書く。ここで文字Cはこの段だけscalar coefficientで、bad vector Cとは別物。
Zはω=3frのtracefree inverseのminusである。
元homogeneous f=a²bを戻しcurlを取ると

`∂c Z=(3E,3H,−2C−G)`。

したがってbad vectorがc方向非零ならE=H零が必要。
log-only sourceのη=1におけるprimitive ψのy²係数、
η=−2におけるr=x²d(ψ/x²)のy²係数はそれぞれH又はEである。
各々ψ=y²,y²+x²（η一）又はψ=y²,y²+x（η負二）の全二型を排除出来る。
source classが各eで存在しても全higherを保った空性であり、標本値を使わない。

この適用ではq(Cbad)とfullw³式からCbad∈C*∂c、κ零を別に確認する必要がある。
R2L γ一源ではqのexact d(b/a)から直ちにその方向になる。
log-only sourceではqの両residue非零は同じ、q=dloga/bの特殊値でも
fullw³式が追加方向を除く（次corollaryで詳述）。

## 6. R3L φ零・r=d(xy)は全空

f=a³、p=q=dloga、r=d(bc/a²)。
raw P0はX=Y零、Z=(0,3b,−3c)、global D2零からK零で、P0自体divergence零。
constant-root character条件はCa零。
canonical strict uv²式は全vectorとしてDZ C零となる。
DZ=diag(0,3,−3)だからCb=Ccも零。
既約quinticにはC非零が必要なので、この全firstconnection sourceは空。

他のr型・他の方向Cへ未確認の結論を拡張しない。
全m三、全G/Eの完了ではない。
