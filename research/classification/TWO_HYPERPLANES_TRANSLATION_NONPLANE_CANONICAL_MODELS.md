# rank二平行移動二全源の構造的canonical評価点

2026-11-08、PRIMARY `E_reduced_nonnormal`。
未計算の全ambient判定に先立ち、二つの完全源それぞれについて一つの構造的な非零parameterモデルを固定する。
別の点を順に選ぶrank farm又は旧normal jetの延長は行わない。

全係数源と全補集合は
`TWO_HYPERPLANES_CUBIC_TRANSLATION_RANK_TWO_NONPLANE_SOURCES.md`。
以下は全sourceの代わりに低次項を零とする仮定ではなく、
その不可約な全源閉包の成分性を最弱の一つの点で判定するための明示模型である。

## 1. L=h、essential四

base座標(h,x)とfiber(u,v,w)を使い、h行及び{x,w}を零とする。

```
{x,u}=x², {x,v}=hx,
{u,v}=3hu−2xv,
{u,w}=−hu+xv−kappa hx,
{v,w}=x²−kappa h².
```

実三次は全kappaで `F=h²u−hxv+x³`、実五次はhxF。
kappa=0のΠ0のPoisson Hamiltonian field V=x∂u+h∂vとmodular fieldは

```
Z=4x∂u+4h∂v−h∂w.
```

両者は可換なPoisson fieldで、上のΠ_kappaは厳密にΠ0+kappa V∧Zである。
したがって全kappaで真のPoisson族であり、形式jetではない。

g_t=diag(1,t,t³,t²,t²)に対し

```
t^(−2) (g_t)_*Π_kappa=Π_(t² kappa).
```

ゆえに複素数体上で全kappa≠0は同じprojective orbit。
唯一の構造的非零代表kappa=1を評価点とする。
元の八parameter全源での値は

```
n4=n7=n9=n10=n15=0, n8=n14=−1, n17=1.
```

全source正常形のfiber変更v_new=v_old−2wにより上の模型へ移る。
元sourceの全八微分もambient方向へ含めるので、この評価点のorbit下界だけで全源を仮定しない。

## 2. L=h−x、essential五

同じcofactor V=x∂u+h∂vを持つ模型を

```
{x,u}=x², {x,v}=hx,
{u,v}=hu−kappa hx,
{u,w}=−hu+xv−xw+kappa x(x−h),
{v,w}=−hu+xv+kappa h(x−h)
```

とする。全kappaで

```
F=(h−x)(hu−xv)+hxw.
```

Π0のmodular fieldはZ=3x∂u+2h∂v+(x−h)∂w。
再びΠ_kappa=Π0+kappa V∧Zであり、V,Zは可換なPoisson fieldである。
g_t=diag(1,1,t,t,t)に対して

```
t^(−1) (g_t)_*Π_kappa=Π_(t kappa).
```

全kappa≠0は同じprojective orbitなので、kappa=1を唯一の非零代表とする。
元の七parameter源の値は

```
n4=n7=n9=0, n8=n14=−1, n10=n15=1.
```

この場合は元のsource座標がそのまま上の模型と一致する。

## 3. 有限全ambient判定の成功条件

二つの元の**全parameter源**を維持して、その全Jacobi・normalizerを再構成する。
各固定点の全150方向D、Euler/GL/scalar及び全source微分Fを取る。
D*F=0、全QQ階数及び非零block、真source切片を保持する。

`rankD+rankF=150` なら、実source次元下界とambient tangent上界が一致し、
不可約な全source閉包が滑らかな全成分となる。
余法線が一又は二なら、同じ点の全Sym² normal→全350/cokerD mapを一回だけ計算する。
そのidealが原点支持なら、真cutとKrull不等式で同じ全成分結論を得る。
実群方向以外のsource方向がある場合、未証明のformal group productを使わない。

非primary支持、余法線三以上、又は実行上限なら、source全体の親又は成分性は未証明として保持する。
高次jet、別kappa代表、別anchor又は部分余核の列へ延長しない。
新成分を数えるには、別作者全ambient確認と全prior・相互の相違監査が別に必要である。
