# rank一short translation全源からlong translation全源への実多項式曲線

2026-11-08、PRIMARY E_multiple_hyperplane。純全source証明。新ambient点、M1、minor、jet又はsolverを使わない。g零／非定数／非零定数と三次共鳴、全Fを同時に保持する。

## 1. 元short源を全て保持する

`MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §5 のshort場はX0=x∂y+∂z。q=(s,w)、r=y−xz、s=xとし、元全sourceは

```
U=u∂s+v∂w, d=div U,
U(g)=d*g, (U+1)H=F*g,
deg g≤1, deg U≤1, deg F≤2, deg H≤3.
```

元の全六括弧は同節(20)であり、gを非零と仮定しない。g零ではH非零がrank四条件、g非零ではrの係数からPfaffianは非零。

## 2. 分母のないactual曲線

λ∈A1に対し

```
s=x−λz, r=y−xz+(λ/2)z²,
Xλ=λ∂x+x∂y+∂z
```

とする。g,u,v,F,H及びそのs,w偏微分を全てこの(s,w)へ代入して、元座標の六括弧を

```
Pλ_xy=−s*u−λ[F+(1+d)r]−r*g_w−H_w,
Pλ_xz=−u, Pλ_xw=λv+g,
Pλ_yz=F+(1+d)r−z*u,
Pλ_yw=xv+zg−r*g_s−H_s,
Pλ_zw=v                                             (1)
```

と定める。全λで元座標の次数は二以下、λの係数も多項式である。λ零は元のshort源に**全六係数で一致**する。

全Jacobiとmodular条件の証明は次の一つのactual多項式自己同型による。補助座標

```
(t,r,s,w)=(z,y−xz+(λ/2)z²,x−λz,w)
```

の逆はx=s+λt,z=t,y=r+st+(λ/2)t²で、Jacobianは全λで−1。これにより(1)は正確に

```
∂t∧Y+Jac(rg+H),
Y=−[F+(1+d)r]∂r+U
```

のpushである。div Y=−1とY(rg+H)=−(rg+H)から全Jacobiが零、modular場はXλ。これは非線形変更から既知PGL族への帰属を推論する操作ではなく、**元座標で二次以下の実Poisson曲線(1)を証明する**操作である。元Pfaffianはrg+Hで非零だから、全λで一般rank四を保ち、射影点が曲線途中で消えることもない。

## 3. λ非零では真のaffine変更でlong源へ入る

λ≠0では次の**affine線形**変更を使う。

```
(tL,Y,q,w)=(x/λ,y/λ,x−λz,w), det=−1/λ.
```

逆はx=λtL,y=λY,z=tL−q/λ。Xλは∂tL+tL∂Yに移り、long正規modular場となる。rL=Y−tL²/2とすれば

```
r=λ*rL+q²/(2λ).
```

long(19)の全parameterは、同じg=g(q,w)を用いて正確に

```
U_L=u∂q+(v+g/λ)∂w, d_L=d+g_w/λ,
H_L=H/λ+q²*g/(2λ²),
F_L=F/λ+[q*u+(1+d)q²/2+H_w]/λ²+q²*g_w/(2λ³).       (2)
```

実際、{q,w}=g、{tL,q}=u、{tL,w}=v+g/λを直接読み戻せる。残る三係数も(1)をこのaffine逆へ代入するとlong(19)と一致する。(2)は全次数制限を保ち、

```
U_L(g)=d_L*g, (U_L+1)H_L=F_L*g
```

を元の二条件から恒等的に満たす。反転したのはλだけで、λ零にこのframeを代入しない。λ零は既に分母なしの(1)で与えられている。

## 4. 全g分岐と閉包への正しい帰結

* g非定数でH=gKならH_L=gK_L、K_L=K/λ+q²/(2λ²)。従って全short非定数g源は `RANK_ONE_LONG_NONCONSTANT_COMPONENT_ENTRANCE.md` の**全long同型源**の閉包へ入る。候補Zの成分性を仮定しなくても、この閉集合への包含は成立する。
* g=g0≠0ならd=0、d_L=0であり、H_Lの三次部はH3/λ。H3零源は全て既受理long定数g/H3零NEW25へ入り、全short同源もNEW25閉包に入る。
* 同じ非零定数gの三次共鳴源は、三次部を消すことなく対応する全long三次共鳴源へ入る。U_Lの線形部分はU1のままなので、全固有値／零係数補集合も保つ。このlong側の帰属が未解決ならshort側もその残余へ戻し、閉じたと数えない。
* g零ではU_L=U,H_L=H/λで、F_Lは依然任意二次式の許容source内にある。従って全short g零源は全long g零源の閉包へ入り、Fを選んだ部分源に切らない。

各場合、λ非零の全元が入る閉集合をCとすれば、実Poisson曲線(1)のλ零もCに入る。これが全short sourceを対応するlong sourceへ送るactual closureの橋である。単なるmodular軌道の退化、元の0/1jetの一致、β像の非優越又は非線形座標同値から帰属を推論していない。
