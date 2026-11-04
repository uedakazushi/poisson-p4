# degree24 の41部分空間に関係式が無いこと：独立係数 minor

2026-11-04。対象は `(F4,p1,p2,p3,p4,N12)` の weighted degree24 monomial
41個だけである。元の beta 点・source 点を評価せず、受理済み七つの完全
多項式 map から係数を取り出す。著者の native integer worker は再実行しない。

独立 engine は有限体 nmod 多項式を使い、変数順を `(B,W,Y,a)`、term order
を LEX に変更する。著者は整数多項式・`(a,Y,W,B)`・DEGLEX だった。各列を
共通分母 `B^12 D^24` へ持ち上げ、保存された41 monomial 行だけの係数を
取り出す。全1681成分を比較し、別の Python elimination で determinant
`23147 mod32003` を得る。全41 monomial の独立な列挙と全columnの整数分母が
p単元であることも確認する。

これは整数係数行列の非零 minor を厳密に示すため、QQ上の列rankも41。
従ってこの41部分空間に NEW source 上で恒等零となる非零な線形結合は無い。
有限体における空集合や source 点の標本から characteristic0 の不存在を
推論したものではない。全係数写像の有限行 minor そのものが証明である。

degree24 の全 Inv24 空間は67次元であり、残り26方向を除外していない。
従って NEW の定義式 degree≥24という既存の下界は**変わらない**。degree24の
不存在、degree26以上、方程式発見は主張しない。高次数への自動 retry は無い。

```sh
. /workspace/cloud-setup/activate.sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_degree24_subspace_review_20261104_v63.py --output work/degree24-subspace-review --work work/degree24-subspace-review-native
```

全体60秒・2GB上限、一回の独立係数 run。コンパイル binary は work に置く。
