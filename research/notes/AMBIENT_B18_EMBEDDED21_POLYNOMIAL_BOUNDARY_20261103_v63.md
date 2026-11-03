# 既存の embedded21 層の多項式境界と B18 原点の21次元形式枝

2026-11-03。W≠0 の全交叉認識の独立受理結果を使う。新しい Poisson 成分や
元の B18 原点の全交叉環を決める結果ではない。

## 1. 全次数の二変数族

```
E(W,λ)=P0−λ h² ∂y∧∂z+W x² ∂z∧∂w,
P0=∂y∧(x²∂h+hx∂z+(h²+z²)∂w).
```

全25発散係数と350生 Jacobi 係数は QQ[W,λ] 上で零。固定 N19 source の
全19係数へ戻して全150係数を検証するので、N への membership に closure
だけを用いる必要もない。W≠0 では rank4、W=0 では rank2となる。

受理済みの群作用をこの全二変数族にも代入すると、厳密に

`r^−1 diag(1,1,r,1,1)_* E(W,λ)=E(W/r,λ)`。

## 2. old embedded prime との実際の係数照合

旧 universal critical-β source の自由係数 k0,k1,k2,k5 の四つを残し、
p=s=r=0、他の k を零にする。potential は

`H=x0*(x1*x4+x3²)+k0*x1²*x3+k1*x1²*x4+k2*x1*x3²+k5*x3³`。

全61正規化座標と、旧 embedded prime の49標的へ戻すと

`v1=−k1, v4=k0, v5=2k2, v7=3k5`

以外は全て零となる。従ってこれは旧 polynomial colon prime の自由四変数
族そのものである。既知の frame の逆を使い、k0=1,k1=λ,k2=0,k5=1/3と
すると、全150物理係数で E(1,λ) に一致する。

W1 の実際の局所 associated prime は、被約 N の極小 prime と一つの
滑らかな21次元 embedded prime である。excellence と忠実平坦完備化により
後者は実際の代数的 associated prime へ降りる。その大域閉包を S21 とする。
GL5×tensor-scalar は連結で D のイデアルを保つので、有限集合 Ass を各々
保ち、S21 を保つ。上の実際の四係数族は局所 embedded prime 上にあり、
恒等式と既約性からその全多項式像が S21 に入る。

よって E(1,λ)、さらに W≠0 の E(W,λ) は全て S21 に入る。S21 は閉なので
W=0 も含み、特に P0∈S21。この議論は新しい大域 prime を別に数えるもの
ではなく、既存の rank4 embedded 層の閉包を P0 へ輸送するものである。

## 3. 原点の21次小行列

P0 の独立な19個の GL 方向と、∂W E、∂λ E を合わせた全150×21行列に、
非零21次小行列がある。全行列と minor を保存する。従って
`GL19局所切断×A²→S21` の像の次元は21であり、この有限逆行列 recipe は
S21 の稠密な21変数有理 parametrization を与える。

P0 でも同じ単位 Jacobian により、この21変数写像は一つの正則な21次元
形式的 image 枝を与える。共通 G19 を除けば embedded support の二次元
形式枝である。しかし S21 の全 germ が滑らか、解析的に既約、あるいは
他の形式枝がないとは結論しない。source の正則な一枝と全 prime の同一視
には、別の接空間上界や proper fiber の検証が必要である。

P_a、Q_b による embedded20 閉包とは次元により異なる。原点の nil support
を曲線だけで説明する模型は、この二次元横断部分を欠く。D 自身の正確な
深さ19/20、nilmodule 全体、全 associated primes の網羅は依然別問題。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_embedded21_polynomial_boundary_20261103_v63.py --output work/B18-embedded21-boundary
```

新しい source jet や Gröbner 基底は実行しない。
