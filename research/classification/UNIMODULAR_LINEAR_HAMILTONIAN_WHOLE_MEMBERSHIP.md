# unimodular最高部の一次Casimir：全低次延長の点別成分被覆

2026-11-08、PRIMARY `E_multiple_hyperplane`。**producer純証明、独立/root査読前。**
Gのprimitive像非優越は既済のまま維持し、ここではEの全点所属を新たに証明する。
新しい点・operator・minor・jet又はCAS計算は実行しない。

## 1. 全source入口と主張

四変数の全二次以下Poisson tensor P=C+L+Qについて

```
div P=0, Pf(P)=k≠0,
Q#dℓ=0 for some nonzero homogeneous linear ℓ
```

とする。sharpは第一添字縮約、div Pのj成分はΣi∂i Pij。
V=P#dℓはaffine-linear、零点を持たず、V(ℓ)=0、divV=0、Lie_VP=0。
`UNIMODULAR_LINEAR_HAMILTONIAN_SUPPORT.md` §§1–3とその独立全45係数査読により、
全(V,ℓ)はconstant、非nilpotent半単純、nilpotent全六対に尽くされる。
cをimageで割った非零類、ℓのchain-covector／他covectorの二軌道を保持する。

**結論。** この全rank四sourceの射影斉次化は、既知有限閉成分の和に入る。
使用する親はCA、T16、T01、rank一long-M1の新25、およびnil22-top A=−1の既22成分。
以下のpolynomial座標はliteral quadratic curveの恒等式を証明するためだけに使い、
非線形同型をPGL所属同値の代わりにしない。

## 2. 非零constant Hamiltonianは全CAへ

Vがconstantならactual affine変更でPxy=a≠0、Pxz=Pxw=0とする。
VのPoisson対称性により全係数はyと独立で、Pf(P)=kからPzw=b≠0はconstant。
divP=0はPyz=A(x,z,w), Pyw=B(x,z,w)に対しAz+Bw=0を与える。
従って全quadratic A/Bを保持したpolynomial primitive H、degH≤3を取り

```
P=a∂x∧∂y+b∂z∧∂w+∂y∧(−H_w∂z+H_z∂w)
```

と書ける。Hのxのみの加法式はtensorを変えない。
元座標のactual quadratic族

```
P_e=(a+e x²)∂x∧∂y+b∂z∧∂w+∂y∧(−H_w∂z+H_z∂w)
```

は全eでPoisson。e非零ではactual x-scale/tensor scaleにより
`MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §4.1の全rawCA源へ入る。
全Hを動かす変換はcubic空間の可逆線形写像で、係数を零へ捨てない。
CAの閉性からe零の元PもCAに入る。
これはconstant-Hamiltonian点がCAのどこかの点に似るという推測ではなく、全六括弧の実族である。

## 3. 非nilpotent affine Hamiltonianは全T16へ

全対の分類により、actual affine frameで

```
V=a(x∂x−y∂y)+∂w,  ℓ=z, a≠0
```

と出来る。semisimple部とconstant部の作用は可換で、一方は半単純、他方は冪零。
Lie_VP=0から両者がPを保つ。従って全係数はwと独立、weight零。
V=P#dzが固定する三括弧と、残三括弧の全weight源は

```
Pxz=−a x, Pyz=a y, Pzw=1,
Pxy=f(z)+c xy, Pxw=x u(z), Pyw=y v(z),
deg f≤2, deg u,v≤1.
```

divP=0からc=0,v=−u、Pf定数からf=f0≠0。
従って全点は二つのcoordinate plane上のsymplectic定数tensor J1/J2と、
それぞれのquadratic Hamiltonian f1=xy、g2=a w+u0 z+u1 z²/2から得るproduct源

```
P=f0 J1+J2+(x∂x−y∂y)∧(−a∂z+(u0+u1z)∂w)
```

に尽くされる。

一般に、別々の二plane上の任意quadratic f,gとそのHamiltonian場Xf,Xgに対し

```
P=A J1+B J2+K Xf∧Xg
```

は全parameterでPoisson、Pf(P)=AB。f/gのnondegenerate Hessian開では
各plane内のactual affine変換でf/gをxy/zwの非零倍＋constantへ移せる。
その開はf/gの全係数空間で稠密であり、全族はirreducibleである。

標準化したproductの全五係数源は

```
{x,y}=p xy+A, {z,w}=q zw+B,
{x,z}=k xz, {x,w}=−k xw,
{y,z}=−k yz, {y,w}=k yw.                         (1)
```

`TWO_QUADRICS_CONFLUENCE_NAMED_MPS_PARENTS.md` の真T16全sourceと、
(A,B,C,D,H)=(x,y,z,w,h)、Euler補正を同資料どおり保って

```
(rho0,rho1,rho2,tau0,tau1)=(q+2k,−p−2k,p−2k,A,B)
```

で全六括弧が一致する。逆は

```
p=(rho2−rho1)/2, q=rho0+(rho1+rho2)/2,
k=−(rho1+rho2)/4, A=tau0, B=tau1.
```

全bindingは定数分母だけで、p/q/k/A/Bの零も含む。従って(1)の全源は閉T16に入る。
一般quadratic f/gのfamilyも上の稠密開とT16の閉性から全T16に入る。
H+相異二quadricの**全pointwise分類は依存にしない**。
named identityを用いない場合も、(1)はA=B=0,p=q=1でSNC点を持つ既約実族なので、
全SNC-meeting component theoremから全MPS有限閉和に入る独立な弱い証明がある。

## 4. 全nilpotent対のための普遍modular化補題

補助polynomial coordinates (t,r,s,l)のJacobianが非零constantで、V=∂t、ℓ=lとする。
V=P#dlから、全Pは

```
P'=−∂t∧∂l+k0∂r∧∂s+∂t∧Y(r,s,l),
k0≠0, div_rsY=0.                                (2)
```

となる。実際、V不変性からt非依存、Pf定数からcomplementary bracket Prsがconstant、
divP'=0からdiv_rsY=0である。Yに任意の高次係数が現れても切り捨てない。

任意のl方向を持たないdivergence零の場D(r,s,l)を取り

```
B'=−l∂t∧∂l+∂t∧D
```

と置く。全三次元Jacobi成分で∂tを二度含む項は零、残るJrsとの交差はdiv_rsD=0で零。
従って任意Yに対しB'自身がPoisson、[P',B']=0、divB'=∂t。
また

```
Pf(P'+eB')=Pf(P')*(1+e l).                       (3)
```

元座標へpushしたBがquadratic以下なら、P_e=P+eBは**元の全sourceに対するactual quadratic
Poisson曲線**で、divP_e=eV、e零が元P、e非零はknown modular sourceに入る。
そのactual射影因子は正確に非零scalar倍 \(h^4(h+e\ell)\) であり、
marked hは全eで重複する。従ってPf(Q_e)=0も成立し、後で適用する全十七rank一source定理や
nil31/22定理の**multiple-H仮定**を満たす。simple-Hの全rank一modular場へ定理を拡大していない。
補助座標そのものから成分帰属は推論しない。

## 5. 五つの非空nil対：全literal Bと全zero補集合

bracket順序は(xy,xz,xw,yz,yw,zw)。第3列は(2)を得る補助座標、最後は元座標のBである。

| V,ℓ | 補助(t,r,s,l) | D | literal B |
|---|---|---|---|
| y∂x+z∂y+∂z, w | (z,y−z²/2,x−yz+z³/3,w) | 0 | (0,0,−yw,0,−zw,−w) |
| y∂x+∂y+w∂z, w | (y,x−y²/2,z−wy,w) | 2(r∂r−s∂s)/3 | ((y²−2x)/3,−2(yz+xw)/3,−yw,−(2z+yw)/3,−w,−w²) |
| y∂x+∂y, w | (y,x−y²/2,z,w) | 0 | (0,0,−yw,0,−w,0) |
| y∂x+∂z, y | (z,x−yz,w,y) | 0 | (−y²,yz,0,y,0,0) |
| y∂x+∂z, w | (z,x−yz,y,w) | 0 | (0,0,−yw,0,0,−w) |

全補助Jacobianは±1である。nil22行では
∂l=∂w+y∂z、∂t=y∂x+∂y+w∂z。
−wV∧∂lのxz成分に生じる−wy²と、Dのpushに生じる+wy²が正確に相殺する。
これにより表のliteral Bは全て元座標で次数二以下であり、source係数に分母を使わない。
各全Pについて(2)–(3)が成り立つので、保存された全45係数kernelの一部だけを選んでいない。

残る六番目、V=y∂x+z∂y+∂w、ℓ=zは、元の全45係数Jacobi/Pf条件の矛盾により空。
`UNIMODULAR_LINEAR_HAMILTONIAN_SUPPORT.md` §3の既独立空性をそのまま保持する。
Vの零点が生じるc∈imageは、P可逆とdℓ非零に反するため存在しない。

### e非零の全既知親

* 第一行：P_e/eのmodular場はnil31-long-topそのもの。実Pfは非零一次式である。
  `MODULAR_NIL31_LONG_SUPPORT.md` の全六源と全Pf式では、このfactor型はA源
  （又は同じA源との交わり）であり、その全sourceはT01に入る。
* 第二行：P_e/eのmodular場はnil22-top。
  全八源のPf分類により非零一次PfはA=−1の全sourceへ入り、その全源は既22成分。
  表のBは同資料(5)のt=0表示にliteral一致する。Bだけのrank二空性を全Pへ外挿していない。
* 後三行：P_e/eのmodular場はrank一long又はshort。全十七long／全shortの点別被覆を使える。
  より具体的には、long座標又はshort座標でg=0、Hは非零affine linear式となりM1へ入る。
  source normalizationで反転するのはeと元Pfの非零constantだけ。
  全e非零点は新long-M1の閉25成分へ入り、shortは既全actual curveで同じM1へ戻る。

各親は閉PGL不変集合であり、P_eは全eで消えない（元Pf k≠0、(3)は恒等零でない）。
従ってe零の**任意の元source**も同じ閉親に入る。
原canonical sourceの自由係数を一つもunitにしておらず、その全零境界を保持する。

## 6. 結合と適用範囲

全affine Hamiltonian対の既 finite分類と§§2–5を合わせ、最高Qに一次Casimirがある
全unimodular・constant nonzero Pf sourceを、有限既知閉和へ点別に被覆した。
R(1,3)の全最高部閉包は、一次Casimirを持つ閉projective incidenceに入るため、
その**全低次延長**にも本定理を適用出来る。generic cubicの正則列を追加仮定しない。

新成分増分は零。全S23/pullback最高部、log M源、因子付き最高部のうちlinear-Casimirでない場合、
零最高部の別直接証明、全A=c=0の完成は、別の正確な全source義務として残す。
元134式/全C6L24Q60との接続はH-Casimir strict liftと全constant-Pf/BV恒等式を用い、
primitive評価のrank不足からこの新E membershipを推測しない。
