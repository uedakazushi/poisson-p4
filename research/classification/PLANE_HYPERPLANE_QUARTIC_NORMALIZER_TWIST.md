# modular Gaの全normalizerから作る実Poisson族

2026-11-07、PRIMARY E。E6 smooth-flex及びA5 smooth-hyperflexの未知親に対する有限命題。
新anchor、高次jet、選んだnormal方向の延長を用いない。

元のfullsourceはh-Casimir代表Πであり、h行は零、全strict Jacobiを満たす。
h=1での実affine modular fieldをNとし、これをh成分零のhomogeneous linear projective fieldへ延長する。
この代表ではNは元Πのdivergenceそのものであり、Poisson vector fieldである。
Eulerを足して作る別のhomogeneous divergence-free liftをNと呼んでいない。
既canonical点ではN²=0、rank N=2で、実全projective stabilizerのGaを生成する。

## 完全なlinear source

全linear vector fields Wについて

\[
[N,W]=aN,\qquad [\Pi,N\wedge W]_{\mathrm{proj}}=0,
\quad a\in\mathbb C
\]

を課す。第一式は25matrix entriesとaのlinear条件であり、第二式も全150係数のJacobi微分に対するlinear条件である。
Euler方向を含む全Wを出発点にし、その射影垂直部分を最後に区別する。
この**完全kernel**をKとする。

各`R_W=N∧W`はPoissonである。実際、[N,W]がNのscalar倍なので
`[N∧W,N∧W]=0`である。
同様に任意のW,V∈Kに対して`[N∧W,N∧V]=0`となる。
各項にはNが二度現れ、又は[N,W],[N,V]が現れるからである。
従って文字通りの全族

\[
\Pi_W=\Pi+N\wedge W,\qquad W\in K
\]

は全parameterでprojective Poissonである。
第一式がstrict normalizer、第二式が全projective Jacobi compatibilityである点を保持する。
各全kernel basisの実matrix・括弧・全350 compatibilityを保存すれば、任意のkernel方向とその積を独立に検証できる。

## 最弱の全成分性への橋

既存の同じcanonical点は全150係数でJacobi derivative rank100、実Euler/GL/scale/source rank48である。
この値と点を変更しない。
全normalizer族の実方向を加えたrankをFとする。

F=50なら、実群作用と上の実族を合わせたmorphismの微分rankが50であり、
raw Poisson像の局所dimensionは50以上となる。
全周囲tangent dimension150−100=50と一致するので、この点の全Poisson local ringはregularであり、
唯一のprojective dimension25 componentを通る。
全族の群飽和閉包は既約かつ同dimensionで、この全componentそのものである。
既受理の全六磁気sourceがこの点の23実軌道閉包に等しいので、その全sourceもこのcomponentに含まれる。
これは単にboundary点が既知族に入ることから逆推論する議論ではない。

F=49ならactual transverse curveが得られた事実とその完全classの一次像までを述べる。
componenthoodやgeneric memberの親をこのrankだけから主張しない。
F=48ならこの全normalizer-twist classの一次像は実群方向内であり、同じ点の部分W probeを続けない。
いずれも未知の全Poisson変形を尽くしたとは主張しない。

計算する場合は両canonical点の一括finite jobで全W kernel、全150/350及び実Fを保持する。
特定のWや新anchorを試す列にはしない。
