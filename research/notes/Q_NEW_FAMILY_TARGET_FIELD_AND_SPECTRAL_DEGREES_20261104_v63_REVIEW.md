# 実第一jet target体の指数1・有理性・spectral次数の独立レビュー

2026-11-04。原著 `Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63`
の全証明と凍結gateを独立に監査した。対象は実15次元第一jet像 `B_NEW⊂W16`。
25次元Poisson成分全体の商体や、ambient W16の商体とは区別する。

## Field inclusionを尽くす

QQ上の実source写像、source derivation、既知の有限τ同値が与える包含は
`E_Q=QQ(B_NEW,Q)^Sp4 ⊂ L_Q=QQ(xi,R,Z0,r,s)`。
これはsourceへの実pullbackであり、任意のhypersurfaceに対しambient invariant fractionsが
全有理不変体を尽くすという誤った一般論を使わない。
独立に受理済みZNのprimitive degree10から
`Kshape_Q=QQ(R,Z0,r,s)=QQ(R,Z0,Z1,Z2)(ZN) ⊂ E_Q`。

全16 source座標によるcocharacter恒等式は任意次数に適用でき、
`F4=xi²*f4`, `J14=xi⁷*h14` となる。f4とh14はshape体に属する。
非零QQ Inv14は既存の3次元Poisson-bracket basisで選べる。
同じ結論は合理的なLie kernelのbase extensionからも従う。
実像の幾何的prime equationのdegree下界26により、そのdegree14 invariantの制限は非零。
従ってf4,h14とも非零で、

```text
xi=(J14/F4³)*f4³/h14 in E_Q.
```

これで全L_Q生成元がE_Qへ入り、`E_Q=L_Q`、指数1、純QQ有理・超越次数5が確定する。
source imageが幾何的に整であることは、QQの座標環を有理source函数体へ入れ、flatなC拡大後も
Cの有理source函数体への単射が保たれることから分かる。同じ包含の証明をC上にも直接適用し、
`C(B_NEW)^Sp4=C(xi,R,Z0,r,s)` を得る。一般の有理不変体がbase changeと可換だと仮定しない。

## 底の違う10・20・40

`K0=C(R,Z0,Z1,Z2)`、`Kshape=C(R,Z0,r,s)` とするとdegreeは10。
xiはshapeから独立なので、`F4=xi²*f4` もshape上transcendentalである。
従って `Kshape(F4)/K0(F4)` はdegree10。
`F4/f4` はF4-adic valuation1を持ち平方でないため、xiを加える次数は2。
よって `E/C(F4,R,Z0,Z1,Z2)` の次数は20。

五つの独立octic係数の体 `P=C(p0,…,p4)` に対し、原著§5の有理逆式を別に記号検査した。
`C(F4,R,Z0,Z1,Z2)=P(F4)` であり、`F4²=400p0+1600p4`。
独立変数のlinear primeのvaluation1からこの拡大も正確に2次なので、`[E:P]=40`。
任意のbase changeとのlinear disjointnessは仮定せず、この二つのvaluationを別々に使う。

projective scalar作用はSp4を除くとxiのweight2 scalingで、shapeを固定する。
その代数的Gm固定体はKshapeなので、projective normalized spectral degreeは10。
QQ上でも同じ代数的characterによる固定体を使い、QQの平方がQQ*を全射するとは言わない。
ZNを加えたlow fieldは `Kshape(F4)` で、全Eはその2次拡大である。

以上はgeneric函数体の次数であり、全base上のfinite/flat写像や全boundary fiberの個数ではない。
元Poisson成分の有理性・全成分尽くし、曲線genus、NEW定義polynomialと正確なdegreeは未解決。
新しい点・prime、Inv14/source26展開、Krylov再実行、巨大なCASは使わない。
