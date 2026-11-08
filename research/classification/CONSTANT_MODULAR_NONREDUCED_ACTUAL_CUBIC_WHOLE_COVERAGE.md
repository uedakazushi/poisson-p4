# constant modular：実cubicが非被約である全sourceのknown被覆

2026-11-08、PRIMARY `E_multiple_hyperplane`、担当 `/root/plane_collision_review`。
**純全source証明。新算術／operator／点／minor／jet／GB無し。**
対象はA=0,c≠0、degf=3の全必要十分源
P=∂t∧Y+Jac(f)、divY=1、Yf=f、degY≤2である。
本稿は**fの実三次斉次化Fが非被約**の場合を全点処理する。
最高項f3の非被約性と混同せず、f3非被約でもF被約の全追加源は残余に保つ。
producerの独立受理は別途必要である。

## 1. 実三次式の全非被約factor

degf=3なので、三次斉次化Fはhで割れない。
非被約三次式はF=L³又はL²M（L,M相異linear）だけである。
L/Mのどちらもhとは異なる。
F=L²Mでaffine因子が平行なら、actual affine座標xで
f=σ x²(αx+β)、αβσ≠0となる。
Yf=fからYx=f/f'=x(αx+β)/(3αx+2β)でなければならない。
分母の根は分子の二根のどちらでもないので多項式にならず、degY≤2に反する。
したがってこの全平行枝は空である。

他の全場合はactual affine座標でf=σx³又はσx²yにできる。
fに現れない第三base座標zをσ⁻¹倍すれば、Jacobian potentialはdetによってσが消え、
κ=1／t座標を保ったままσ=1へ移せる。この可逆尺度はσ≠0でのみ用い、
σ=0はdegf=3という入口の外である。

## 2. cube f=x³：全quadratic Yの別marked H

Yf=fからYx=x/3が全域で従う。Yy,Yzは任意の二次以下のまま保持する。
元h-Casimir斉次liftのdivergenceはN_h=−h∂t。
canonical divergence-free liftはΠbar=P−E∧N_h/5である。
column-sharp規約でP#dx/x=(h/3)∂t、したがって

```
N_x=5 Πbar#dx/x=(2/3)h∂t.
```

xはactual D=h²x³の三重hyperplaneなので、x-Casimir liftの全源を取れる。
x=1の新affine chartではN_xはrank一nilpotent linear field、定数部零。
全Yの係数を特別化せず、`MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §4の
rank一/c-image全source定理に入る。
従って全cube sourceはCA/CB/既nil31-DDHの有限閉和へ入る。
この結論にY affine又はgeneric cubic Hamiltonianという仮定は要らない。

## 3. double+simple f=x²y：全十係数の必要十分表示

Yf=fの全式は2yYx+xYy=xy。
mod xとgcd(x,y)=1からYx=xA1、Yy=y(1−2A1)、A1は全affine式である。

```
A1=λ0+λx*x+λy*y+λz*z.
```

divY=1を代入すると∂zYz=λ0+3λy*y+λz*zであり、全解は

```
Yx=xA1,
Yy=y(1−2A1),
Yz=λ0*z+3λy*y*z+λz*z²/2+K2(x,y),
```

K2は任意の二次式六係数となる。四λとK2のA10全体で全条件は恒等的に成立する。
逆も上の全係数比較により成立し、追加Jacobi compatibilityを捨てていない。
全パラメータでsigned Pf=f=x²y、actual D=h²x²yでrank四である。

## 4. 新x∞の全modular型とcritical零補集合

base順序(x,y,z)のJac(f)から、column-sharpで

```
P#dx/x=A1_h∂t+x∂z,
A1_h=λ0*h+λx*x+λy*y+λz*z.
```

N_h=−h∂tとcanonical gaugeを用いると、x-Casimir liftのdivergenceは

```
N_x=[(5λ0−1)h+5λx*x+5λy*y+5λz*z]∂t+5x∂z.
```

x=1のaffine chartでlinear部Aはt出力の一行だけ、
定数部は5λx∂t+5∂zである。
A≠0ならrankA=1、A²=0、定数部はimage(A)=C∂tの外にある。
λz≠0ならA c≠0でlong型、λz=0ならshort型。
同じ全rank一/c-nonimage源へactual affine Jordan正規化出来るので、
`RANK_ONE_TRANSLATION_WHOLE_COMPONENT_COVERAGE.md` の全十七long源／全short曲線被覆を使う。
全F、全period、全zero境界を持つ有限既知閉和Kへ入る。
ここで新fixed点のfactorやβ像上界を所属証明に使わない。

残るA=0は正確にλ0=1/5、λy=λz=0である。
これは全A10必要十分sourceの閉部分であり、λxとK2は一切制限しない。
そのsourceでλ0だけを変えるactual polynomial Poisson族はA≠0開へ入り、
Pf=x²yは不変。Kは有限閉和なので、criticalの全λx,K2もKへ戻る。
自分自身のA0 branchへ循環的に所属を仮定していない。

## 5. 正確な全scope

F=L³、非平行L²M、平行L²M空の全場合と全係数境界を戻した。
よってA0,c≠0、degf=3、**actual F非被約**の全源はknown閉成分で覆われる。
actual affine/PGL/frameと非零tensor尺度を戻しても同じ結論が保たれる。
新成分増分零、smooth／local unique／全分類は非主張。

まだ未処理なのは、最高f3がdouble-line+line又はtriple-lineなのに
actual Fが被約となるsourceである。
既に処理したf3被約全源、degf≤2全源、及びA=c=0のunimodularh⁵源とは区別する。
