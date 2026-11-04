# Full content の独立最大性とsource shape拡大次数10

2026-11-04。原著 `Q_NEW_FAMILY_TWO_SHAPE_CONTENT_GCD_20261104_v63` の
70項content C、全15 cofactor、891項Gを独立に検証する。
元のraw determinant identityと、sのgeneric map degrees3/5は別の独立packetで受理済み。

## 全係数の積と異なる最大性certificate

独立したPython整数係数演算で、全15群 `Q_ij=C*q_ij` と literal全体 `Q=C*G` を検査する。
さらに既知linear factorを反復乗算し、47,114全係数について

```text
raw=25*2^68*r^22*Z0*(r−1)^16*T^15*C*G
```

を独立に確認する。原著のnative multivariate gcdやresultantは再実行しない。

最初の係数群(0,3),(5,0)はZ0を含まない非零多項式なので、全群の任意の共通因子は
`QQ(R)[r,Z0]` でZ0次数0である。C除去後の全群をZ0係数へ分け、
**QQ(R)[r]上の別のEuclid/extended Bézout**で選択した有限係数の恒等式

```text
Σ U_i(r,R)*q_i(r,R)=1
```

を作る。全U_iの分子・分母と、使った係数labelを保存し、恒等式自体を再確認する。
R-only分母はgeneric coefficient fieldの単元であり、共通分母を掛ければ非零D(R)を右辺に
持つ多項式certificateとなる。Z0の係数を取ることは、Z0=0の空性を一般へ持ち上げる
操作ではない。共通因子が既にZ0-freeと分かっているから正当化される。
これで全contentの最大性を別方式で証明する。Rの特別値の幾何へは拡張しない。

## Primitive implicitizationとexact degree

`k=QQ(r,R,Z0)` 上でZ1(s),Z2(s)のmap degreesは3,5。
従って `k(Z1,Z2)=k(s)` であり、sによるplane curve parametrizationはbirational。
個別分子・分母が互いに素という既受理gateによりbase pointはなく、primitive resultantは
その既約曲線のimplicit equationを一乗で与える。
今回、全(Z1,Z2)係数contentを除いたので、Gauss lemmaから
`G∈QQ(R)[r,Z0,Z1,Z2]` は既約である。

受理済みoctic/spectral baseの優越性と
`p4/F4²=1/1600−Z0/4`, `p3/F4²=Z0(R+1)−1/400`
により `R,Z0,Z1,Z2` は独立である。Gはrを実際に含み、そのr次数は10。
もう一度Gaussを用いれば `G∈QQ(R,Z0,Z1,Z2)[r]` はrの最小多項式。
sはrとこの底から有理的に回収できるので、**既知sourceのshape体について正確に**

```text
[ QQ(R,Z0,r,s) : QQ(R,Z0,Z1,Z2) ] = 10.
```

これは未調査boundaryを全て捨てるglobal/properness結果ではない。
N12/F4³の原始性、full NEW target不変体とsourceの指数1、NEW定義方程式・最小次数を
証明したものでもない。ξの回収に必要な全source homogeneityは別gateである。
新点、prime、resultant、PRS、source26、原著multivariate gcdの再実行は使わない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_two_shape_content_review_20261104_v63.py --output work/v62-review/two-shape-content-independent
```
