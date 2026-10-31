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

## v27：既知三成分・全S61局所検証・一次図式照合

2026-10-31。T39/T40もW16へ優越的で、一般に有限であることを証明した。
既知三成分の指定点は全S61内で滑らかで、一次射影はétale。
十分一般の一次部分には少なくとも三つの相異なる孤立被約な延長がある。
全一般ファイバーの尽くしや有限性、次数は未証明。

```sh
.venv/bin/python code/elliptic_weight1_dominance.py --output reports/new-elliptic.json
.venv/bin/python code/s61_known_fiber_audit.py --output reports/new-known-fiber.json
.venv/bin/python code/verify_s61_fiber_certificate.py --output reports/new-fiber-independent.json
.venv/bin/python code/toric_coframe_fiber.py --output reports/new-coframe.json
.venv/bin/python code/verify_toric_coframe_independent.py --output reports/new-coframe-independent.json
.venv/bin/python code/mps_catalog_audit.py --output reports/new-mps-crosswalk.json
```

- `elliptic_weight1_dominance.py`：一次資料のq(5,1)/q(5,2)一般有理族、Jacobi恒等式、枠の正規化、16×16微分小行列式。原文(21)の誤った負号を使わず(18)とJacobiで確認した正号を用いる。
- `s61_known_fiber_audit.py`：元の134式を直接微分。T01/T39/T40の61座標、復元、Jacobi、全階数45、固定ファイバー階数45を確認。
- `verify_s61_fiber_certificate.py`：解析的な微分コードを使わず、重み2で二次以下という事実と厳密な中心差分で6,075個の保存行列成分を検証。
- `toric_coframe_fiber.py`：中心化した五超平面をR∈GL4で復元。16変数方程式と自由S5商でトーリック開部分のファイバーを記述。境界・他成分は含めない。
- `verify_toric_coframe_independent.py`：逆射影変換の直接微分と一般六パラメータで全双ベクトルを比較し、S5生成元を検証。
- `mps_catalog_audit.py`：Table 2の全40図式をS5同値込みで一次資料照合。`--source-html`で元HTMLからの再抽出も検証する。オフライン用の事実抽出は `mps_table2_diagrams.json`。

全51階数四成分のS61座標表示が完成したという意味ではない。
楕円型の全一般パラメータ61座標の保存と、三射影の一般次数は残る。
