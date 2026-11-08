# 定数modular source Iの全点をM11全源閉包へ送る実曲線

2026-11-08、PRIMARY E、producer plane_collision_review。純手証明。
**全f=zw／source Iは既存のrank一long・M11 whole sourceのPGL飽和閉包Z_M11に入る。**
Z_M11のambient成分性／既知親はなお未解決なので、本結果を既知成分へのK被覆とは呼ばない。
新成分増分は零。新operator、point、rank、kernel、profile、GB、jetは実行しない。

## 1. 完全A10源のactual normalizer

κ=1、f=zw、Jac orientation(u,z,w)で

```
Yu=u*(-βz+γw)+K2(z,w),
Yz=z*(αu+βz+γw+δ),
Yw=w*(1-αu-βz-γw-δ)
```

が全sourceである。K2は全六係数。
α≠0でactual affine shear

```
u'=u+(βz+γw+δ)/α, z'=z,w'=w,t'=t
```

を用いる。det一、inverseは同じlinear/constant shiftを引く。
全括弧を変換するとnormal mixed partsはYz=αu'z、Yw=(1-αu')w、
Yu'はu'非依存の全quadratic Kbarとなり

```
Kbar=K2+(β²/α)z²-(γ²/α)w²
          +(βδ/α)z+[γ(1-δ)/α]w.
```

βγδの零補集合を初めから落としていない。
Kbarのz,w一次係数をl_z,l_wとすると、actual linear t shear

```
t'=t+l_z*z-l_w*w
```

はこの二項を消す。逆はminus。Jac(f)のHamiltonianによりYuが
`Yu-l_z*z-l_w*w`へ変わり、Yz,Ywは変わらない。

残るKbarをk0+p z²+q zw+r w²と書く。
αpr≠0のopenで

```
u'=a*u,t'=a*t,z'=b*z,w'=c*w,
tensor scale=1/a,
a=α,b²=αp,c²=αr
```

とする。全scale非零、C上の有限平方根を選ぶだけでactual GL/tensor scaleの逆を持つ。
結果は

```
P_tu=F=n+z²+mzw+w²,
P_tz=u*z, P_tw=(1-u)*w,
P_uz=z, P_uw=-w, P_zw=0,
n=αk0, m=αq/(b*c).
```

m,nに追加unit条件は無い。generic sourceの二moduliを無断で固定していない。

## 2. M11の全raw sourceからの実曲線

`RANK_ONE_TRANSLATION_WHOLE_FINITE_EXHAUSTION.md` のM11を使う。
source座標は(x,y,z,w)、H=τzw、U=λz∂z+μw∂w、μ=-1-λ、τ≠0。
全quadratic F_old(z,w)を残した元六括弧は

```
P_xy=-F_old,
P_xz=λz, P_xw=μw,
P_yz=z*(λx+τ), P_yw=w*(μx-τ), P_zw=0.
```

これはr=y-x²/2を非線形PGLと見なした表示ではない。
d=λ+μ=-1なので、元fullsourceのPxy=-F_old-(1+d)rのr項が正確に消えた式である。
全λを含み、λ零／μ零も有効な元whole sourceである。

任意の同じF=n+z²+mzw+w²に対し、ε≠0で

```
λ=1/ε, μ=-1-1/ε, τ=ε^-3,
F_old=F/ε³,
t=ε²*y, u=ε*x+1/ε, z=z,w=w
```

とする。mapはactual affine GL4、homogeneous GL5はhを固定しdet=-ε³。
逆はx=(u-1/ε)/ε、y=t/ε²。tensor scaleは一、Euler gaugeは零。

Leibnizとinverseを全六括弧へ代入すると正確に

```
P_tu=F,
P_tz=u*z,
P_tw=[1-(1+ε)u]*w,
P_uz=z,
P_uw=-(1+ε)w,
P_zw=0.
```

例えばPtz=ε²z(λx+τ)=z(εx+1/ε)=uz。
Ptw=ε²w(μx-τ)=w[-(1+ε)(u-1/ε)-1/ε]=[1-(1+ε)u]w。
Ptu=ε³P_yx=Fである。他三式も表示通りで、負冪は一つも残らない。
全εでPfaffianはzw、全h-Casimir liftのactualDはh³zwである。
元sourceではPf=-τzwであり、det=-ε³との積も同じ恒等式を与える。

元M11 sourceの全Poisson性と実GL変換によりε≠0でPoisson。
表示は全εでquadraticかつε多項式なので、元Jacobiの閉性からε=0もPoisson。
ε=0は§1のcanonical source Iそのものである。

## 3. 全sourceへの閉性と未解決点

Z_M11は全M11 sourceのPGL飽和像のZariski閉包であり、閉PGL不変である。
したがって全ε≠0の曲線とそのendpointがZ_M11に入る。
§1のactual normalizerにより、αpr≠0の全source Iも同じZ_M11へ入る。
fixedfの元sourceはA10で、tensor mapは全parameter上regular／非零。
Z_M11の閉逆像がdense openを含むので、α=0,p=0,r=0を含む**全A10**へ結論が延びる。
単位normalizerのformulaへ零を代入せず、元の全sourceの閉性を使う。
全元κ、base affine frame、modular cとHのPGL frameを戻してもZ_M11は保たれる。

この包含はEの二義務を同じwhole親問題へまとめるものである。
M11の旧exact25像上界／一点D98から成分性又は既知成分所属を推論しない。
source IIは別の全NEW25 theorem、source IIIはなお独立に未解決である。
本稿は新count、全constant-modular E完成又は全分類完成を主張しない。
