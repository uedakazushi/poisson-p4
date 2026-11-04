# Actual CBの四つの正則g商：原点値はすべて厳密に0

2026-11-04。全186相対入力に残っていた四つの未評価原点値を、一つの既存unit57系の
三次二変数jetによって決定した。新しいmembership、colon、nilradical、全Dのdepthは
計算していない。既存の凍結DAGの`None`値を変更するpacketでもない。

## 1. 既受理の入力と、今回求めるもの

`AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63`および独立REVIEWは、actual CBの
同じphysical base R0において、次の九つが原点を保つ正則座標であることを保証する。

```
(z1,...,z7,r,s)
 =(f18,f19,f28,f29,f60,f67,f72,j46,j95).
```

除外したzero odd generatorはf63である。z=0への制限を
`g(r,s)=f63(X(z=0,r,s))`とする。四つのzero-weight graph row
`alpha=64,68,73,74`のactual odd coefficientを同じ逆座標で制限したものをbbar_alphaと書く。
既受理のzero-image algebraによるmembership証明は、全次数で

```
bbar_alpha = g*h_alpha,    h_alpha in QQ[[r,s]] algebraic
```

を与える。ここでh_alphaは既存full DAGの四つの`g`除算node
`330164,330289,330524,330594`であり、当時は原点値のみを未評価として保存した。
今回の対象は四つの数`h_alpha(0)`だけである。

actual conductor multiplicity4、mu(I_CB)=10と既存の九座標からgは非零である。
同じphysical wにおけるconductor quotientは

```
QQ[[w,r,s]]/(g,wr,ws)
   = QQ[[w]] fiber_product_QQ QQ[[r,s]]/(g).
```

従って受理済みのG19 productを除いたmultiplicity4は`1+ord(g)`であり、ord(g)=3となる。
この理由は旧analytic modelの座標を新しいphysical wと同一視することに依存しない。
今回の三次jetもord(g)=3を直接確認する。

## 2. 一つのpointed57系だけの計算

元の41315-node source DAGから、CBの48本のsource方程式、上の九座標、gと四つの
bbarに必要な9978 ancestorだけを使う。未知数は元のsource coefficient48個と
physical zero base9個の計57個である。方程式は

```
source equations48=0,
F9=(0,...,0,r,s).
```

原点の有理57×57 Jacobian Jは

```
det J = -3456/625 != 0.
```

JとJ逆行列は証明書に全保存し、両者の積が単位行列であることを照合した。
`QQ[r,s]/(r,s)^4`において、同じ固定J逆行列を使って次数1、2、3の順に
同次Newton補正を求めた。各次数の残差の最低次項に対する線形化は原点Jであり、
これによってpointed inverseの三次jetが一意に決まる。
最後に元の57本の方程式すべてが全10 monomialで成立することを照合した。
g除算node自体やfull339669-node DAG、16個のinverse cloneは評価していない。

## 3. 保存した全係数と結論

monomial順序を

```
1,r,s,r^2,rs,s^2,r^3,r^2s,rs^2,s^3
```

とするとgの三次jetは

```
g = -24r^3-48r^2s+18rs^2+54s^3 + O((r,s)^4)
  = -6(r-s)(2r+3s)^2 + O((r,s)^4).
```

四つのbbarは、この順序の全10係数がすべて0である。
既知の正則恒等式`bbar=g*h`の最低次成分は
`bbar_3=g_3*h(0)`である。g_3のlex順で最初の非零係数r^3は-24なので、
四つともh(0)=0と一意に決まる。選択した係数だけでなく全四つの三次係数も
`bbar_3=g_3*0`と一致することを確認した。

| physical graph row | full DAG g商node | 正確な原点値 |
|---:|---:|---:|
|64|330164|0|
|68|330289|0|
|73|330524|0|
|74|330594|0|

この推論で、jetから全次数membershipを推測してはいない。正則性と全次数の除算恒等式は
前packetの入力であり、有限jetから読み取ったのは商の定数項だけである。
商が恒等的に0であること、gがこの三次多項式に等しいこと、必要な高次jetの上界は主張しない。

## 保存と計算範囲

- Source: `code/ambient_B18_CB_g_division_origins_20261104_v63.py`。
- Evidence: `reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63/degree3-certificate.json`。
- 同folderのreportに全source hash、実測時間、メモリを保存。

唯一の実計算は0.708555秒、最大RSS131972KiB、証明書96626bytes。
Pythonの既存Fraction演算を使い、インストールや環境変更はしなかった。
実計算前の起動時にpython-flint不在とmetadata keyの取り違えを検出・修正したが、
いずれも方程式評価前に停止し、jetの再試行・次数増加は行っていない。
GB、新しい点、source/core再生成、p逆元、ungauged R9系、高次jetは使わない。
全Dのdepth、T3のcolonやsocleを決定するには、別の実際のmodule計算が必要である。
