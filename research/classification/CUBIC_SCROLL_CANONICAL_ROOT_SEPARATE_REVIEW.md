# smooth cubic scroll canonical q2-root：全空性の独立検証

2026-11-07、別作者 `/root/quartic_review`。**PASS（純証明）**。
producer `CUBIC_SCROLL_CANONICAL_ROOT_EXCLUSION.md`、SHA256
`819efc1448fe9fa39c5d5f999c9c4bc4cffff31e799eaafdea315a577678bbca` を全文読解した。
新CAS/source/minor/Jacobi/GB/jet計算はない。

smooth S(1,2)⊂P4、actualD被約既約quintic、S⊂SingD、restriction generic rank2、
normal quadratic Gにprimitive canonical root ℓ0（dq2のsaturation）がある全源を対象とする。
normal rank1/2、全可約/非被約R、root collision、codim2 basepoints、vertical成分を許す。
**この全源は空**。directrix必要条件の残五pair中ℓ0を含む四組とrank1 A0を排除し、
A=2C0+5fのrank1/rank2二型及びcubic coneは残す。

## canonical liftの三恒等式

canonical divergence-zero homogeneous quadratic lift ΠはGL5-equivariantに一意に取れる。
projective Jacobi E∧[Π,Π]=0とBV divergenceを組み合わせ、Euler weight0/ambient dimension5から
strict [Π,Π]=0。任意raw代表をstrict/div0と仮定しない。

α=ι_(Π²)Ωはclosed、係数degree4なのでd(ιEα)=5α。
actual quintic DはιEαの非零定数倍、rank4線形代数Π#α=0よりΠ#dD=0。
q=q2=x2x4−x3²、X=X_qとすると

```
X(D)=0, X(q)=0, div X=0.
```

Hamiltonian divergence零はΠのdiv零及びHessian対称性からも直接従う。
以後Euler gaugeを変更してこの三恒等式を壊さない。

## rank3 quadricの全finite cover

Q=(q=0)のaffine coneはA=C[x0,x1,z,w]のμ2 invariant ring、
x2=z²,x3=zw,x4=w²、involution(z,w)→(−z,−w)。
weights=(2,2,1,1)、T=x0w−x1zはweighted degree3のprime。
S逆像idealは(zT,wT)なので、height1 supportは正確にT。
vertex側のembedded境界もT内に含む。genericSではcover étaleである。

dqがnormal rootなのでD|Qはgeneric Sでmultiplicity≥3。
primitive saturationがC0でdqを割っても、genericSでは同じroot lineでありこの含意は変わらない。
UFD AでT³がDbarを割る。q∤DよりDbar≠0。
weighted degDbar10/degT³9から**全許容形**は

```
Dbar=T³ell, ell=az+bw, (a,b)≠(0,0).
```

一係数零も保持する。multiplicity≥4ならweightでDbar0、即q|Dとなり除外される。

## intrinsic root-line preservationとnonreduced Hartogs

parallel quadratic G/σはgeneric root lineを保ち、rank1でも唯一root lineを保つ。
chart x0≠0で正しいintrinsic numeratorは

```
Hi=x0{q,xi}−xi{q,x0}−2q{x0,xi} ∈ (q)+I_S².
```

raw {q,xi}だけへ置き換えない。
cover上でx0X(xi)−xiX(x0)がT²倍。
gcd(x0,T)=1から、全indicesでxiX(xj)−xjX(xi)がT²倍となる。
これはgeneric identityをprime Tのdivisibilityで全域へ延長したもので、chart/vertexを除いたままにしない。

B=A/(T²)はdim3 hypersurfaceなのでCM、originでdepth3。
五つのxi=(x0,x1,z²,zw,w²)の共通零はoriginだけで、各xiはBでnonzerodivisor。
上のminors identityよりpunctured Spec Bの各D(xi)上でX(xi)/xiが一つのregular functionにglueする。
H0_origin(B)=H1_origin(B)=0のHartogsからB全体の元Λへ延長。
normalityではなく**nonreduced CM depth**だけを使っている。
weightΛ4<weightT²6なのでAのweight4 polynomialに一意にliftし、

```
X(xi)=Λxi+T²ci, ci constant.
```

Xq0はc2w²−2c3zw+c4z²=0を強制し、c2=c3=c4=0。
X(z)=Λz/2、X(w)=Λw/2がgeneric étale liftを与え、branch locusも越えてregular。
したがって全cover derivationは

```
X=Λ Etilde+T²V,
Etilde=x0∂x0+x1∂x1+(z/2)∂z+(w/2)∂w,
V=c0∂x0+c1∂x1.
```

頂点以外のconstant方向を無証明で捨てず、Xq0から消した。

## Casimirと正しいvolume divergenceの3対2

Etilde(Dbar)=5Dbar、Vell0なのでXDbar0は
`5Λ+3T V(T)=0`。
一方Qのcanonical residue volumeはcover上2 dx0∧dx1∧dz∧dw。
実際dx2∧dx3/x2=2dz∧dwで、branchに新しいzero/poleを足さない。
ambient divX0/Xq0からこのvolumeも保存される。
divEtilde3、EtildeΛ2Λ、divV0なので
`5Λ+2T V(T)=0`。
両式の差がTVT0、integral AでVT=c0w−c1z0、従ってV0及びΛ0。
この係数差はnormalizationの一般性やsmooth conductorを仮定したものではない。

よってX|Q0、全{q,xi}はqで割れ、Qはprojective Poisson hypersurface。
rank4 symplectic locusのsmooth点にPoisson hypersurfaceは存在しないからQ⊂Dが必要。
D被約既約degree5に反する。

全root rankと全R/vertical/codim2境界を保つ空性であり、
旧kernel16/弱いG28等の未閉数値結果を成功へ読み替えていない。
残A25二型、cone、全G/Eの尽くしはこのreceiptから結論しない。
