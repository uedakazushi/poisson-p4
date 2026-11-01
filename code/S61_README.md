# S61 の現行座標・一次部分・既知族

2026-11-01 / v36-s61。原本は `legacy/v02/` に置いたまま、以下の独立実装と
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
v29時点の関係はN_k=50 e_k d_kである。v30でN_k=3000、d_k=1、e_k=60と決定した。

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
t−1となり、この制限上のファイバーは一つの被約点である。五変数全体の一般長さは後続v30で3000と決定した。

## v30: 普遍的切断、交叉、既知五点の三平面

`q51_universal_section_20261101.py` と `q52_universal_section_20261101.py` は全16変数上の45二次式、全134代入恒等式、直接Jacobi検算を保存する。
同名の `_review_20261101.py` は別の復元・線形代数で点検する。
`elliptic_section_intersection_20261101.py` は旗円錐と原点冪零部分、同名reviewは全Hilbert級数も独立確認する。
`toric_universal_three_plane_20261101.py` とreviewは普遍的中心と既知五点の三平面を構成し、`s61_five_cover_product_20261101.py` は積代数と冪等元を保存する。
**三平面への全ファイバー包含は未証明。** `curve_invariant_smooth_counterexample_20261101.py` は誤った輸送保存則への反例であり、追加方程式として使用しない。
全スクリプトには `--output work/新しい名前.json` を指定し、保存済み証明書を上書きしない。
詳しい適用範囲と独立点検は `research/notes/README.md` と各JSONの入力ハッシュを参照。

## v31: 二捩れによるトーリック成分の記述

`toric_jacobian_two_torsion_20261101.py` は優越的六次元断面上の全結合則恒等式を確認する。
`toric_e2_scalar_inverse_20261101.py` は相互逆座標とtraceを独立に検証。
`primitive_curve_discriminant_20261101.py` は8次jet行列と24次平方自由直線制限を、
`toric_c9_polynomial_20261101.py` は階数低下を避ける射影直線のBézout証明書を保存する。
`toric_polynomial_e2_coordinate_20261101.py` はXをβ,qの多項式で復元し、実際のT01のΔ非零上の全域同型を与える。
`toric_missing_boundary_20261101.py` と `toric_boundary_intersections_20261101.py` は欠落する8次元軌道と境界交叉を確認する。
このT01の記述を全S61の尽くしと混同しない。

## v32: 無限遠と特殊ファイバー

- R0無限遠の二平面型と非二元型: `S61_R0_INFINITY_20261101`。対応する code / reports / research/notes に再現資料を保存。
- 二つの楕円五次曲線の直交双対性: `NULL_CORRELATION_DUALITY_20261101`。対応する code / reports / research/notes に再現資料を保存。
- 一般固定因子頂点ファイバーの根基: `S61_R0_VERTEX_RADICAL_20261101`。対応する code / reports / research/notes に再現資料を保存。
- 非二元無限遠点の完成局所環: `S61_NONBINARY_FORMAL_RING_20261101`。対応する code / reports / research/notes に再現資料を保存。
- 旗型無限遠点の完成局所環: `S61_FLAG_FORMAL_RING_20261101`。対応する code / reports / research/notes に再現資料を保存。
- 実際のT01零一次部分の次元10〜12: `TORIC_ZERO_JET_FIBER_20261101`。対応する code / reports / research/notes に再現資料を保存。
- T02一次像の定義式次数は26以上: `MPS_T02_INVARIANT_OBSTRUCTION_20261101`。対応する code / reports / research/notes に再現資料を保存。

## v33: 零一次部分の正確な次元と一般無限遠の帰着

- 滑らかな全領域での既知五延長の分離: `KNOWN_FIVE_SMOOTH_LOCUS_20261101`。対応する code / reports / research/notes に証明・再現資料を保存。
- 全頂点パラメータの一様な根基: `S61_R0_UNIFORM_VERTEX_20261101`。対応する code / reports / research/notes に証明・再現資料を保存。
- 一般一次部分のR0被約無限遠の尽くし: `S61_R0_GENERIC_EXHAUSTION_20261101`。対応する code / reports / research/notes に証明・再現資料を保存。
- 一般βの零モジュラー無限遠をR0へ帰着: `S61_C0_INFINITY_20261101`。対応する code / reports / research/notes に証明・再現資料を保存。
- 非平方二平面型の局所次元3: `S61_BINARY_LOCAL_DIMENSION_20261101`。対応する code / reports / research/notes に証明・再現資料を保存。
- 実際のT01零一次部分の次元は10: `TORIC_ZERO_JET_DIMENSION_20261101`。対応する code / reports / research/notes に証明・再現資料を保存。
- 滑らかな標準曲線を持つ非零モジュラー無限遠の反例: `S61_MODULAR_INFINITY_COUNTEREXAMPLE_20261101`。対応する code / reports / research/notes に証明・再現資料を保存。

- ABは零一次部分の唯一の最大次元成分: `TORIC_ZERO_JET_TOP_COMPONENT_20261101`、独立点検 `TORIC_ZERO_JET_TOP_COMPONENT_REVIEW_20261101`。
- トーリック境界のstrict安定化群障害: `TORIC_STRICT_STABILIZER_20261101`、独立点検 `STABILIZER_OBSTRUCTIONS_REVIEW_20261101`。
- 滑らかな全領域の有限安定化群: `SMOOTH_CURVE_FINITE_STABILIZER_20261101`、独立点検 `STABILIZER_OBSTRUCTIONS_REVIEW_20261101`。

## v34: 零一次部分の尽くしと無限遠の局所環

- 零一次部分のR0全点はABに入る: `TORIC_ZERO_JET_R0_20261101`。対応する証明・コード・証明書と独立点検を参照。
- 最小型の全零一次部分をABへ帰着: `TORIC_ZERO_JET_MINIMAL_20261101`。対応する証明・コード・証明書と独立点検を参照。
- 非零トレースの全点をABへ帰着: `TORIC_ZERO_JET_NONZERO_TRACE_20261101`。対応する証明・コード・証明書と独立点検を参照。
- 実際のT01零一次部分の被約台はAB: `TORIC_ZERO_JET_EXHAUSTION_20261101`。対応する証明・コード・証明書と独立点検を参照。
- ABの正規proper incidenceと全ファイバー: `AB_INCIDENCE_20261101`。対応する証明・コード・証明書と独立点検を参照。
- ABの特異点集合を完全決定: `AB_INCIDENCE_20261101`。対応する証明・コード・証明書と独立点検を参照。
- 二平面型の横断長さ65と非零の斉次化変数: `S61_BINARY_TRANSVERSE_ALGEBRA_20261101`。対応する証明・コード・証明書と独立点検を参照。
- 全135二次斉次化のSp4同変性: `S61_HOMOGENIZED_EQUIVARIANCE_20261101`。対応する証明・コード・証明書と独立点検を参照。
- 一般βに残る二つの最小冪零型無限遠曲線: `S61_MINIMAL_NILPOTENT_INFINITY_FAMILY_20261101`。対応する証明・コード・証明書と独立点検を参照。
- 最小型の平方零完成局所環と飽和による除去: `S61_MINIMAL_INFINITY_FORMAL_RINGS_20261101`。対応する証明・コード・証明書と独立点検を参照。

## v35: 一般無限遠の全型帰着と局所飽和

- ABの原点外の正規性と二つの接錐: `AB_ONE_VARIABLE_LOCAL_20261101`。証明・証明書・独立点検を参照。
- 被約連結ファイバーだけから正規性は従わない: `REDUCED_FIBERS_NORMALITY_COUNTEREXAMPLE_20261101`。証明・証明書・独立点検を参照。
- 冪零(2,2)型の一般無限遠からの除外: `S61_NIL22_INFINITY_20261101`。証明・証明書・独立点検を参照。
- 正則冪零型と最小型開層の除外: `S61_OTHER_NILPOTENT_INFINITY_EXCLUSION_20261101`。証明・証明書・独立点検を参照。
- 全非冪零モジュラー型の除外: `S61_NONNILPOTENT_INFINITY_EXCLUSION_20261101`。証明・証明書・独立点検を参照。
- 最小型無限遠の全点に対する頂点条件: `S61_MINIMAL_INFINITY_VERTEX_NECESSITY_20261101`。証明・証明書・独立点検を参照。
- 一般β上の最小型被約無限遠の尽くし: `S61_MINIMAL_VERTEX_EXHAUSTION_20261101`。証明・証明書・独立点検を参照。
- 一般一次項上の被約無限遠を全型で尽くす: `S61_GENERIC_REDUCED_INFINITY_EXHAUSTION_20261101`。証明・証明書・独立点検を参照。
- 最小型M0のaxis成分の一般ファイバーの除去: `S61_MINIMAL_AXIS_INFINITY_20261101`。証明・証明書・独立点検を参照。
- 二平面型の一様な二次・三次障害と長さ下界65: `S61_BINARY_UNIFORM_RANK_20261101`。証明・証明書・独立点検を参照。
- 非二元型の二つの共鳴交点と飽和除去: `S61_NONBINARY_RESONANCE_NODES_20261101`。証明・証明書・独立点検を参照。

## v36: 二平面型の自由性、非平坦な交点と新しい例外曲線

- 二平面型の自由階数65と非零平方零のt: `S61_BINARY_FORMAL_RING_20261101`。明記した証明・係数記録・独立点検を参照。
- 最小型M1の三乗零の尾: `S61_MINIMAL_CUBIC_TAIL_20261101`。明記した証明・係数記録・独立点検を参照。
- 四追加直線の平方零完成局所環: `S61_MINIMAL_BASELINE_THICKENINGS_20261101`。明記した証明・係数記録・独立点検を参照。
- 非二元u49=0点を通る二つの最小型枝: `S61_MINIMAL_ZERO49_BRANCHES_20261101`。明記した証明・係数記録・独立点検を参照。
- 非平坦な三枝合流点の一般長さ7と飽和除去: `S61_NONBINARY_ZERO49_ENDPOINT_20261101`。明記した証明・係数記録・独立点検を参照。
- AB正規化の次数二に新生成元はない: `AB_ORIGIN_DEGREE_TWO_20261101`。明記した証明・係数記録・独立点検を参照。
- 一般βにも残る二平面型の階数降下曲線: `S61_BINARY_RANK_DROP_FAMILY_20261101`。明記した証明・係数記録・独立点検を参照。
- 滑らかな標準曲線でも二平面型階数は下がる: `S61_BINARY_RANK_DROP_SMOOTH_CURVE_20261101`。明記した証明・係数記録・独立点検を参照。
