# dependent normalnil・trace零の全超平面還元：別作者監査

2026-11-07。reviewer `exceptional_coverage`、producer `generic_exhaustion`。
入力 `PLANE_RANK_ZERO_DEPENDENT_TRACE_ZERO_HYPERPLANE_REDUCTION.md`
SHA256 `c8ef8dadb1f3a14cdee166e714093364dc41c7ade6cfd101b93684c51809523f`。
**全 necessary coefficient・actual hyperplane・matrix 還元を独立に手計算して PASS。**
新CAS・点計算・全残余のG閉鎖を主張しない。

入口は可換conormalを持つrank零planeのcanonical35/30/9 source、global V!=0、
normalnil、trace field T=0。全30/9 higherを保つ。
Π0=zV wedge partial_w+alpha*z² partial_z wedge partial_w、V=A*x,trA=0,A!=0。
Π1の一項w²C wedge partial_zに着目する。

## 1. 全係数中の forcing を直接再導出した

bilinear strict JacobiをΠ0とこのΠ1項から計算すると、base indices i,jについて

```
J_ijz=2zw(V_i*C_j-V_j*C_i),
J_ijw=-w²(V_i*C_j-V_j*C_i),
J_izw=zw²*(AC+2alpha*C)_i.
```

最後の AC はCがVを微分する項、2alpha*CはCがalpha*z²を微分する項であり、
どちらも同じ符号である。Euler wedgeのbase i,j,z,w成分は

```
x_i J_jzw-x_j J_izw+z J_ijw-w J_ijz
=zw² * [x wedge(AC+2alpha*C)-3(Ax) wedge C]_ij.
```

従ってproducerの係数、特に−3と2alphaの符号が一致する。
zL,wM、zwB、zwE、z²Dの同quadruple成分はz³又はz²wのnormal monomialであり、
zw²には寄与しない。Π2及びhigher括弧は別weightなのでこの式を変更しない。
他係数を零とするsliceを選んだ必要条件ではない。

## 2. C零の補集合は actual D factor で除ける

全canonical sourceで contraction with dzはz倍の項と符号を除きw²Cだけである。
Π2はbase-base、vertical normalはalpha*z²なので別のz=0項は無い。
C=0なら全π#dzがz倍であり、z=0は実際のprojective Poisson hypersurfaceとなる。
そのtangent dimensionは三、skew rankは二以下であるから、実Pfaffian五次にz因子が生じる。
rank4かつ実Dが被約既約quinticというGの仮定に反する。
global divergence liftやformal normal coneを用いた推論ではない。

## 3. C非零の全matrix恒等式と残余

forcingへx=Cを代入すると4C wedge AC=0。
AC=mu*Cと書き直すと、((mu+2alpha)I−3A)xの全imageはC-lineに入る。
あるlinear form ellにより

```
A=((mu+2alpha)/3)I+C tensor ell,
ell(C)=-(mu+2alpha)
```

となる。二つ目はtrA=0による。
さらにAC=mu*Cを戻すとmu=−2(mu+2alpha)/3、従って

```
mu=-4alpha/5,
A=(2alpha/5)I+C tensor ell,
ell(C)=-6alpha/5.
```

alpha!=0ではC-lineとquotientの固有値が−4alpha/5,2alpha/5で異なる。
rank-one補正のminimal polynomialもこの相異なる二根を割るので、Aはsemisimple double型であり、
mixed型のJordan補集合は残らない。
alpha=0ではA=C tensor ell,ell(C)=0、A!=0からrank一nilpotentのみ。

受理済み全五Jordan表へ戻すと、regular semisimple、mixed、regular nilpotentは全alphaでG空。
double semisimple H=diag(1,1,−2)ではquotient eigenvalue1よりalpha=5/2だけ。
Cは−2のdistinguished lineに入る。
rank一nil Nではalpha=0だけ、CはimageNの非零vectorである。
Vのscale/Jordan正規化は元のnoteのactualGL3/normalGL2/scalarと逆変換へ明示的に依存する。

残るのは全H・K=(5/2)z²と全N・K=0の二source。
前者のsymbol base上界六、後者五も、旧trace零表のalpha固定による正しい値である。
これらの全30/9 higher解・G非優越・E帰属は本必要条件だけでは結論しない。
T非零、normalss、scalar V零、又は全高重複度の尽くしへ拡張しない。
