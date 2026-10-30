# S61 の現行座標・一次部分・既知族

2026-10-30 / v26-s61。原本は `legacy/v02/` に置いたまま、以下の独立実装と
証明書から研究を再開する。完全分類・素分解・新しい成分を主張する実装ではない。

## 座標規約

`s61_model.py` は5×5反対称行列の斉次二次双ベクトルを入力とする。
150係数の順序は `(i,j)` の辞書順、次に二次単項式の
`combinations_with_replacement(range(5),2)` 順である。
Euler垂直成分24方向を有理数体上のRREFで除き、v02の保存された61自由係数を使う。
`parameters[k]` は旧データの `u(k+1)`。重み1と重み2は連続したブロックではない。

`normalize(P)` の前提は固定点 `[1:0:0:0:0]` での値が
`J = d1∧d2+d3∧d4` に一致すること。別の点を選ぶ操作とこの値の正規化は呼出側で行う。
4個のperspective係数を有限の射影変換で消す。
`reconstruct(u)` は保存済み55復元式から150係数を戻す。
`verify_point` は復元の一致・134式・独立なprojective Jacobiを確認する。
非斉次・非交代行列は拒否する。浮動小数点を研究証明に使わず、有理数または有理関数を渡す。

`verify_equations` は一般パラメータでの分母を共通化し、分母を払った多項式を
QQ上の疎多項式環で比較する。分母の零点へは結論を延長しない。

## 再現

setup後、新しい出力先を指定して実行する。

```sh
.venv/bin/python code/verify_s61_adapter.py --output reports/new-adapter.json
.venv/bin/python code/verify_s61_families.py --output reports/new-families
.venv/bin/python code/s61_weight1_audit.py --output reports/new-weight1-gauge.json
.venv/bin/python code/weight1_structure.py --output reports/new-weight1-structure.json
.venv/bin/python code/toric_weight1_dominance.py --output reports/new-toric-dominance.json
.venv/bin/python code/weight1_rank_one_fiber.py --output reports/new-rank-one-fiber.json
```

各CLIは出力上書きを避ける。S61の再生成コマンド、原本との比較、ハッシュは
`reports/S61_REBUILD_20261030.json` を参照。

追加の再現：

```sh
.venv/bin/python code/ambient_weight1_image.py --output reports/new-ambient-image.json
```

CA/CBの全一次部分の共通6次元線形族、任意の非退化定数項の正規化、
P3上の階数6 incidence、像の次元9を厳密に検証する。

## 一般パラメータ族

トーリック族は6個のlog Poisson係数 `a,b,c,d,e,f` をすべて残す。
これは旧CSVのbiresidue行列とは異なる。chartは `a*(a*f-b*e+c*d) != 0`。
`toric_family` の後 `normalize` を適用すると、この開集合上のS61写像が得られる。
出力 `toric_T01_map.json` は61個の有理関数を保存する。
単一の6パラメータ写像を16次元S61成分の全パラメータ表示とは呼ばず、残存群作用を加える。

もう一つの族はv02のalmost-abelian二変数族で、特殊化を含む全 `alpha,beta` について
恒等式を検証する。この族には今回、MPS40内の個別成分IDを割り当てていない。

## 一次部分と残余問題

保存済みgaugeと、自然なprimitive gauge `tau(L)=L12+L34=0` は異なる。
有限perspective変換で移り合うが、重み2変数も変わるため、一次変数だけを交換しない。
`s61_weight1_audit.py` はこの有限変換を独立に検算する。
`weight1_structure.py` はprimitive閉一次2形式の16次元空間、群作用、階数1軌道を扱う。
`toric_weight1_dominance.py` は一般トーリック族の一次部分を独立に構成し、
6パラメータ微分とSp4作用を合わせた階数16の厳密証明書を作る。
`weight1_rank_one_fiber.py` は既存61変数から階数1の固定fibreを取り出す。
fibreの接空間次元38は実際の次元・成分数・被約性を意味しない。
