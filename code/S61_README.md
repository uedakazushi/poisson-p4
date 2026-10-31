# S61 の現行座標・一次部分・既知族

2026-10-31 / v29-s61。原本は `legacy/v02/` に置いたまま、以下の独立実装と
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
全一般ファイバーの尽くしや有限性、次数は未証明。v28では下記の通りT01の次数を3と決定した。

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
楕円型の全一般パラメータ61座標の保存と、二つの楕円射影の一般次数は残る。

## v28：トーリック次数3と全模型の無限遠

```sh
.venv/bin/python code/toric_algebra_fiber.py --output reports/new-toric-algebra.json
.venv/bin/python code/toric_algebra_fiber_review.py --output reports/new-algebra-review.json
.venv/bin/python code/toric_algebra_degree.py --output reports/new-toric-degree.json
.venv/bin/python code/toric_degree_review.py --output reports/new-degree-review.json
.venv/bin/python code/s61_infinity_flags.py --output reports/new-infinity-flags.json
.venv/bin/python code/s61_infinity_flags_review.py --output reports/new-infinity-review.json
```

- `toric_algebra_fiber.py`：固定一次部分をもつtrace-free可換積を20個の対称三階テンソルで表示。96二次式で結合性を課し、非零トレース判別式上でGL4/S5と同定する。Poisson双ベクトルとノルム因子を復元。
- `toric_algebra_fiber_review.py`：一般40係数から線形写像を再構成し、上ブロックを含む64恒等式、トレース形式と文字による復元を独立に確認。
- `toric_algebra_degree.py`：96式から普遍的一次帰結18本を保存。二変数化した全イデアルの商が指定点で長さ3の被約代数となり、三点とも半単純であることを確認。二次先頭係数の可逆性が開集合上の有限性を保証し、有限étale近傍により一般次数3を証明。
- `toric_degree_review.py`：別の密行列核計算から96式全体を消去し、全証明書の行操作・商代数・Bézout恒等式も検算。
- `s61_infinity_flags.py`：未飽和二次斉次化の旗incidence、55復元式の重み2線形性、Bianchi20関係の係数次数を確認。最終証明書は `S61_INFINITY_FLAGS_20261031_final.json`。これを真の射影閉包やアフィン正次元性の証明としない。

`s61_infinity_flags_review.py`は元のモデル実装をimportせず、20個のBianchi恒等式、
55復元式、方向双ベクトルの直接Schouten計算、全134先頭式を独立確認する。

全一般S61には少なくとも五つの孤立被約な延長がある。一般次数3はT01の射影に限り、
楕円二成分・未知成分・例外層の尽くしを含まない。

## v29：楕円型五変数模型と全ファイバーの消去

```sh
.venv/bin/python code/elliptic_six_variable_fiber.py --output reports/new-elliptic-five.json
.venv/bin/python code/elliptic_five_variable_review.py --output reports/new-elliptic-five-review.json
.venv/bin/python code/s61_fiber_quadratic_span.py --output reports/new-quadratic-span.json
.venv/bin/python code/verify_s61_macaulay_degree3.py --seconds 150 --entry-limit 50000000 --output reports/new-macaulay3.json
```

楕円型のコード名は中間の六変数還元に由来する。最終関数
`fixed_jet_five_variable_equations` はλ,p1,…,p4の五変数について30等式を返す。
共変frameの行列式とI4,I14・正規化分母を可逆化し、Sp4行列と倍率を一意に復元する。
分子を巨大な多項式に展開せず、有限個の加減乗除からなる処方を保存する。
一般解数N_kは未計算で、成分の一般次数d_kとの関係はN_k=50 e_k d_k、e_kも未決である。

全ファイバーについては134個の二次先頭部の非零小行列式をQQで保存した。
三次Macaulayの完全計算はF32003に限られ、6164×17296行列の階数6073、
追加二次式一つ、一次帰結なしを確認。C++コンパイラが必要で、一時フォルダへビルドする。
時間・成分数上限に達した実行は部分結果として記録する。
標数零での追加式、一般次数、全尽くしの証明に読み替えない。

全Macaulayの観測から普遍的帰結への持ち上げは、次の証明書で完了した。

```sh
.venv/bin/python code/s61_weight4_syzygies.py --output reports/new-universal-weight6.json
.venv/bin/python code/s61_weight6_review.py --output reports/new-weight6-review.json
.venv/bin/python code/elliptic_fiber_line_probe.py --output reports/new-elliptic-line.json
```

70本の重み4先頭式の一次関係は36次元で、35次元のSchouten--Bianchi部分の外の
関係からG6を作る。全61変数で9257項、重み6、q次数は高々2。元の134式の
多項式係数結合なのでイデアルを変えず、指定点で二次先頭階数135を持つ。
独立点検は有理核の上界と有限体階数の下界を組み合わせてQQの階数を証明し、
全G6と135次小行列式を確認する。指定flag代表上の先頭制限は恒等的に零である。

楕円型の直線制限はλ=2,p=(t,2,3,4)のみを扱う。QQ(t)上の二つの必要式のgcdが
t−1となり、この制限上のファイバーは一つの被約点である。五変数全体の解数は未計算。
