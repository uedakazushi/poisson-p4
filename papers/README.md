# 公表用原稿の構成案と証明資料

**原稿・PDF・研究台帳は v61-s61、2026-11-03。固定二座標の有限射影候補、境界の剰余次数、非平方 a=0・d=0 層の一般除去を統合。一般五点の尽くしは未完了。**

[一次ジェット還元の日本語原稿](s61_first_jet_reduction_ja.tex)（[PDF・110ページ](s61_first_jet_reduction_ja.pdf)）は、
primitive gauge、T01/T39/T40の優越性と全S61のétale点、楕円型五変数模型、T01の一般次数3、
少なくとも五つの一般延長、トーリック代数模型と無限遠の旗、階数1固定fibre、
CA/CBの共通像、全模型の普遍的な重み6帰結を扱う。MPS40の全図式と内部IDも一次資料へ照合済み。
標準曲線の判別式、Jacobianの非零二捩れとの同型、実際のT01のΔ非零上の有限エタール性、全域の飽和イデアルと零ファイバーABを扱う。
全一般アフィンファイバーの次元上界2と元二座標への有限射を使い、任意の固定値で一度の t² コロンから零次元射影候補を得る。未飽和の三次元基本サイクルの全次数2300を、アフィン長とは区別する。指定平方旗の二因子除去と c=2a の全パラメータ多項式模型、非平方 K=a=0・K=d=0 の全一般β層の混合標数による除去も記す。この二層で非平方 marked 境界を尽くしたとは主張しない。
最終原稿・PDFの検証範囲と再ビルドは [最終同期記録](../reports/MANUSCRIPT_SYNC_20261103_FINAL.md)。[初回同期記録](../reports/MANUSCRIPT_SYNC_20261103.md)の109ページ版は対応するソース・PDFとともに保存した。
内部研究草稿であり、全成分の尽くし、新規性の網羅的調査と投稿準備は未完。
再現手順は [現行コード](../code/S61_README.md)。

ここは公表用原稿を再構成するための入口である。現時点での本文は、`sources/v02.zip` から `sources/v24.zip` に保存した研究草稿にある。セットアップ後は `legacy/v02/` から `legacy/v24/` に展開される。下の三本は投稿原稿の**候補構成**であり、完成原稿、査読済み論文、完全分類の達成を意味しない。各版の証明・証明書は保存してあるが、全依存関係の独立査読と新規性の網羅的調査は未完了である。

研究全体の現在地はルートの [README](../README.md)、評価資料は [研究ノートの入口](../research/notes/README.md)、引用の適用範囲は [REFERENCES.md](../research/REFERENCES.md) を参照する。機械可読の原稿索引は [SOURCE_INDEX.json](SOURCE_INDEX.json) にある。

## 共通の対象と主張の範囲

基礎体は複素数体とする。主対象は

$$
\mathcal P=\{[\pi]\in\mathbf P H^0(\mathbf P^4,\wedge^2T_{\mathbf P^4}):[\pi,\pi]=0\}
$$

であり、非零定数倍で射影化するが、$\operatorname{PGL}_5$ 商は取らない。以下の成分次元はすべてこの空間での次元である。「成分の尽くし」「全射影軌道の分類」「非被約スキーム構造の決定」を区別する。

原稿で現在主張されているのは、rank $\le2$ の被約 locus の8成分、全 $\mathcal P$ の一般 rank2 成分6個、および全成分数の下界

$$
40+6+11=57
$$

である。一般 rank4 成分は少なくとも51個であり、その尽くしは未解決である。57を正確な総数として書かない。

## 三本の候補

| 候補 | 主題 | 主原稿 | 現在の状態 |
| --- | --- | --- | --- |
| A | rank2 Poisson 構造の被約成分の尽くし | v14 の安定分岐と全8成分定理 | 証明付き研究草稿。2026-10-30 に v14 の10検証群を再現。独立査読・引用監査を経て統合する |
| B | SNC 領域と交わらない追加成分と一般点での非被約構造 | v3・v5・v6・v7・v8 | 11成分の構成・識別・局所環計算が分散。今回の評価では旧版の全大規模証明書を再実行していない |
| C | 正規な実際の退化因子をもつ構造の分類 | v14 および v16–v22 | 正規因子なら $q_{5,1}$ 型という閉点分類を主張。長い入力鎖の統合と独立監査が必要 |

分割は固定しない。独立した定理と証明を読者が追えること、既知研究との差が明確であることを基準に調整する。A・B・Cの完成を、全 rank4 分類が終わるまで待つ必要はない。

## 候補 A：rank2 の成分の尽くし

### 取り込む結果

| 結果 | 原資料 |
| --- | --- |
| 安定な飽和接層をもつ分岐は実際に rational Jacobian 型となる | [v14 stable_rank2.tex](../legacy/v14/sections/stable_rank2.tex)、`thm:stable-rank-two-actual-rational` |
| 被約 rank $\le2$ locus の8成分と、一般 rank2 の全成分6個 | [v14 rank2_global_eight.tex](../legacy/v14/sections/rank2_global_eight.tex)、`thm:rank2-eight` |
| 非安定分岐・因子的零点をもつ分岐の被覆、二つの境界成分の rank4 への変形 | 同上、および同節が指定する保存済み証明・計算 |
| 今回の再現結果と監査範囲 | [global_audit.txt](../audits/2026-10-30/global_audit.txt)、[v14 verification_summary.json](../audits/2026-10-30/v14_recheck/verification_summary.json) |

8成分は $\operatorname{Rat}(1,1,3)$、$\operatorname{Rat}(1,2,2)$、$\operatorname{Log}(1,1,1,2)$、LPB、$\operatorname{Aff}(5)$、$\operatorname{Aff}(4,1)$、toric rank2 $\operatorname{Log}(1,1,1,1,1)$、$\operatorname{Aff}(3,2)$。最後の二つは一般 rank4 成分の境界に入る。

### 統合時の作業

1. 定義、canonical lift、飽和接層の安定性、原始的な閉3形式の対応を一つの設定へ揃える。
2. v14 の安定分岐を中心に据え、reflexive sheaf と locally free sheaf、元の5切断と全切断空間を区別した証明を保持する。
3. 非安定分岐に用いる Corrêa–Muniz の本文 Theorem 4.1 の仮定を明記する。文献の stable branch の主張を、v14 の rational 型の尽くしと同一視しない。
4. 8族が被覆すること、8族が成分であること、うち6個が全 Poisson 空間の成分であることを別々に証明する。
5. 独立な検証器と証明書を対応づけ、どの幾何学的論点が機械検算の外にあるかを記す。

全軌道、境界の全交叉、rank2 locus の冪零構造はこの論文の主定理に含めない。解決済みの安定 rank2 分岐を、旧版の未解決リストを理由に再開しない。

## 候補 B：追加11成分と非 SNC の幾何

下表は既存草稿で全成分と主張されている11個の索引である。名前は内部識別子であり、文献固有の分類番号ではない。

| 内部識別子 | 射影次元 | 主証明の所在地 |
| --- | ---: | --- |
| smooth-double-quadric | 22 | [v3 主原稿](../legacy/v03/p4_poisson_completeness_v3.tex)、smooth double quadric の成分定理 |
| triple-hyperplane | 22 | [v5 triple_hyperplane.tex](../legacy/v05/triple_hyperplane.tex) |
| $\mathcal C_A$ | 28 | [v6 ambient_families.tex](../legacy/v06/ambient_families.tex)、`thm:AB` |
| $\mathcal C_B$ | 28 | 同上、`thm:AB` |
| nil22-open | 23 | [v6 nilpotent22.tex](../legacy/v06/nilpotent22.tex)、`subsec:quarticcomponent` |
| zero-log | 23 | [v6 zero_log_component.tex](../legacy/v06/zero_log_component.tex)、`thm:zero_log` |
| M | 23 | [v7 minimal_thickness.tex](../legacy/v07/minimal_thickness.tex)、`thm:minimal-full` |
| A0 | 23 | [v7 nil22_v7.tex](../legacy/v07/nil22_v7.tex) |
| D± | 23 | [v7 dplus_fullp4.tex](../legacy/v07/dplus_fullp4.tex)、`thm:dplus-full`。符号を別成分として数えない |
| C+ | 24 | [v8 cplus_component.tex](../legacy/v08/sections/cplus_component.tex) |
| C− | 23 | [v8 cminus_component.tex](../legacy/v08/sections/cminus_component.tex) |

総数の整理は [v8 count_and_frontier.tex](../legacy/v08/sections/count_and_frontier.tex)、2026-10-30 の評価は [global_audit.txt](../audits/2026-10-30/global_audit.txt) にある。

### 統合時の作業

1. 11成分について、一般点の表示、実際の退化因子、全接空間、実際の族の次元、成分性の証明、相互識別の不変量を同じ形式の表へ揃える。
2. MPS の40成分との区別は、一般点の因子が SNC でないことだけで済ませない。保存原稿で用いる閉じた incidence 条件が成分全体で SNC locus との交わりを排除することを再確認する。
3. 一般点での非被約構造は、接次元の超過だけで結論しない。高次障害、実際の jet、完成局所環の同定を原稿から取り込む。
4. 被約な退化因子をもつ追加成分が、MPS の問題提起にどう応えるかを一次文献と照合する。現在の文献全体での新規性確認は未完了なので、「初めての反例」などの優先権の主張を先に書かない。
5. v3–v8 の全証明書を再監査し、各成分の条件と例外値を本文だけでも追える形にする。

v23 の SS と NN は同じ既知成分 $T_{37}$ に同定されている。v24 の $Q\in T_{37}\cap\mathcal C_A\cap\mathcal C_B$ や26次元共通族は、追加成分を増やす結果として数えない。これらの境界解析を候補Bの後半に入れるか独立ノートにするかは、本文の分量と主定理の見通しから決める。

## 候補 C：正規退化五次因子の分類

### 主張と入力鎖

[v22 setup.tex](../legacy/v22/sections/setup.tex) の `thm:normal` は、一般 rank4 の Poisson 構造について、**実際の因子** $D=V(\pi\wedge\pi)$ が正規なら、滑らかな楕円正規五次曲線 $C$ が存在し、$D=\operatorname{Sec}(C)$、$\pi$ は $q_{5,1}$ 型であると主張する。これは複素閉点の分類である。$D_{\mathrm{red}}$ の正規性だけを仮定する結果でも、全 Jacobi スキームの被約性の結果でもない。

| 段階 | 主資料 |
| --- | --- |
| 正規五次因子からの一意復元、Gorenstein Jacobian 条件、determinantal なスキーム表示 | [v14 normal_quintics.tex](../legacy/v14/sections/normal_quintics.tex) |
| 滑らかな inactive 成分の制約 | [v16 smooth_curves.tex](../legacy/v16/sections/smooth_curves.tex) |
| inactive isotropy の冪零性と一次法線項の帰着 | [v17 nilpotent_isotropy.tex](../legacy/v17/sections/nilpotent_isotropy.tex)、[v18 setup.tex](../legacy/v18/sections/setup.tex) と同版の各制約節 |
| inactive 直線の排除 | [v19 setup.tex](../legacy/v19/sections/setup.tex) と三重・四重型の各節 |
| 特異 inactive 成分の排除、全 $\widetilde E_6$ 型の完結、$\widetilde E_8$ 四次支持への帰着 | [v20 secant.tex](../legacy/v20/sections/secant.tex)、[inactive.tex](../legacy/v20/sections/inactive.tex)、[e6.tex](../legacy/v20/sections/e6.tex)、[e8.tex](../legacy/v20/sections/e8.tex) |
| 超平面内の支持の排除と八配置への帰着 | [v21 setup.tex](../legacy/v21/sections/setup.tex)、[space_geometry.tex](../legacy/v21/sections/space_geometry.tex)、[trees.tex](../legacy/v21/sections/trees.tex) と関連排除節 |
| 八配置の全排除 | [v22 主原稿](../legacy/v22/p4_poisson_completeness_v22.tex)、[verification.tex](../legacy/v22/sections/verification.tex) |
| 独立監査の到達点 | [geometry_strategy.txt](../audits/2026-10-30/geometry_strategy.txt) |

### 独立監査の優先箇所

- v17–v19：特異支持全体が一本の直線である場合だけに成立する旧結果を、可約支持中の一成分へ使っていないか。一定の法線旗、三重・四重型の尽くし、例外枝の排除を重点的に確認する。
- v20–v22：被約支持の特異点・接触・多重点を含む網羅性と、必要 first-jet 条件が全候補へ適用できること。多項式恒等式の再現だけではこの部分の証明にならない。
- 文献：Pym の purely elliptic 分類を、正規性だけからの分類と読み替えていないこと。今回の追加議論が外す仮定を Introduction で明記する。

v22 末尾にある「一般階数二が残る」という記述は、v14 の到達点を反映していない。旧原稿を保存したまま、統合稿の研究状況では修正する。

## 全分類の研究と、この三本の関係

全一般 rank4 成分に対応する $S_{61}$ は [v2 主原稿](../legacy/v02/p4_poisson_completeness_v2.tex) に構成されている。16個の重み1変数と45個の重み2変数、134方程式からなる実際の代数モデルである。45変数の二次断面や v24 の端点の局所枝だけで、この全モデルの成分を尽くしたことにはならない。

候補A・B・Cは、全分類を目指す研究を支える独立した成果候補である。三本を書き終えたことを完全分類の終了条件にはしない。一方、全分類が未完であることを理由に、既存の証明を公表可能な形へ整える作業を先送りしない。

## 原稿を取り込む際の規則

1. `legacy/` は原資料として扱い、訂正は新しい原稿へ書く。版ごとの ZIP と照合できる状態を保つ。
2. 引用する定理ごとに、元の版・ファイル・定理ラベル・仮定・検証記録を `SOURCE_INDEX.json` で追えるようにする。
3. 書き直した日本語原稿は LuaLaTeX を用いる。記入欄や書込み用ボックスを追加しない。
4. 一つの候補の全本文を出力する前に、主定理と必要な補題の依存図、引用監査、独立検算の不足を確定する。
5. 未確認の部分は未確認のまま表示する。原稿の記述を強めるには、証明または検証記録の追加を必要とする。

v35では一般被約無限遠の全型帰着、ABの原点外の正規性、二平面型の三次障害上界、非二元共鳴交点を統合した。全分類は未完了。

## v36: 二平面型の自由性、非平坦な交点と新しい例外曲線

- 二平面型の自由階数65と非零平方零のt: `S61_BINARY_FORMAL_RING_20261101`。明記した証明・係数記録・独立点検を参照。
- 最小型M1の三乗零の尾: `S61_MINIMAL_CUBIC_TAIL_20261101`。明記した証明・係数記録・独立点検を参照。
- 四追加直線の平方零完成局所環: `S61_MINIMAL_BASELINE_THICKENINGS_20261101`。明記した証明・係数記録・独立点検を参照。
- 非二元u49=0点を通る二つの最小型枝: `S61_MINIMAL_ZERO49_BRANCHES_20261101`。明記した証明・係数記録・独立点検を参照。
- 非平坦な三枝合流点の一般長さ7と飽和除去: `S61_NONBINARY_ZERO49_ENDPOINT_20261101`。明記した証明・係数記録・独立点検を参照。
- AB正規化の次数二に新生成元はない: `AB_ORIGIN_DEGREE_TWO_20261101`。明記した証明・係数記録・独立点検を参照。
- 一般βにも残る二平面型の階数降下曲線: `S61_BINARY_RANK_DROP_FAMILY_20261101`。明記した証明・係数記録・独立点検を参照。
- 滑らかな標準曲線でも二平面型階数は下がる: `S61_BINARY_RANK_DROP_SMOOTH_CURVE_20261101`。明記した証明・係数記録・独立点検を参照。

## v37: 普遍非二元型の除去と接続層の一般化

- 二平面型の階数降下を越える自由階数65: `S61_BINARY_RANK_DROP_LOCAL_20261101`。係数証明書と独立点検を参照。
- 二平面型の二つの端点層を全一般βで除去: `S61_BINARY_ENDPOINTS_20261101`。係数証明書と独立点検を参照。
- 最小型八接続点の一般階数とtの冪零性: `S61_MINIMAL_JUNCTION_BRANCHES_20261101`。係数証明書と独立点検を参照。
- 最小型八接続層の全一般βファイバーを除去: `S61_MINIMAL_JUNCTION_GENERIC_STRATA_20261101`。係数証明書と独立点検を参照。
- AB正規化の次数三にも新生成元はない: `AB_ORIGIN_DEGREE_THREE_20261101`。係数証明書と独立点検を参照。
- 非二元三共鳴層の一次像を除外: `S61_NONBINARY_RESONANCE_EXCLUSIONS_20261101`。係数証明書と独立点検を参照。
- 非二元xzw軌道の全一般β無限遠を除去: `S61_NONBINARY_CANONICAL_EXHAUSTION_20261101`。係数証明書と独立点検を参照。
- 三つ目の二平面型端点と非一般の交差: `S61_BINARY_U4_ENDPOINT_20261101`。特殊標本の追加共鳴と一般層を区別する。
- 最小型M0/M1の全四孤立軸表示を除去: `S61_MINIMAL_ALL_AXIS_REMOVAL_20261101`。元134式の全階数45と相対正則性を用いる。

## v38: 純二重根軌道の全除去と二平面型の全階数降下表

- 二平面型の十八次元階数降下族: `S61_BINARY_HYPERSURFACE_RANK_DROP_20261101`。明記した依存・係数証明書・レビュー範囲を参照。
- 非平方二平面型の全階数降下表: `S61_BINARY_RANK_DROP_COVER_20261101`。明記した依存・係数証明書・レビュー範囲を参照。
- 二つの十七次元二平面型層の自由階数65: `S61_BINARY_CODIM_TWO_FORMAL_20261101`。明記した依存・係数証明書・レビュー範囲を参照。
- 非対角主曲線公式とM1階数一頂点の除去: `S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101`。明記した依存・係数証明書・レビュー範囲を参照。
- M0階数一主曲線とA1の接続点: `S61_MINIMAL_RANKONE_JUNCTION_20261101`。明記した依存・係数証明書・レビュー範囲を参照。
- M0階数一A1の全有限部分を除去: `S61_MINIMAL_RANKONE_A1_REMOVAL_20261101`。明記した依存・係数証明書・レビュー範囲を参照。
- 純二重根型の普遍平方零局所環: `S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101`。明記した依存・係数証明書・レビュー範囲を参照。
- 純二重根型の見掛け三因子を除去: `S61_NONBINARY_DOUBLE_ROOT_ARTIFACTS_20261101`。明記した依存・係数証明書・レビュー範囲を参照。
- 二重根と二最小尾の一般階数五接続: `S61_NONBINARY_DOUBLE_ROOT_ENDPOINT_20261101`。明記した依存・係数証明書・レビュー範囲を参照。
- 二重根型の四非平坦層の除去: `S61_NONBINARY_DOUBLE_ROOT_TORSION_DIVISORS_20261101`。明記した依存・係数証明書・レビュー範囲を参照。
- 二重根と混合尾の自由階数五の接続: `S61_NONBINARY_DOUBLE_ROOT_MIXED_JUNCTIONS_20261101`。明記した依存・係数証明書・レビュー範囲を参照。
- 純二重根二表示の全一般β無限遠を除去: `S61_NONBINARY_DOUBLE_ROOT_EXHAUSTION_20261101`。明記した依存・係数証明書・レビュー範囲を参照。
- AB次数四は `AB_ORIGIN_DEGREE_FOUR_EXPLORATION_20261101` に未完の探索として保存。

## v39: 四混合曲線軌道と最小型の有限基線・半単純尾の除去

- R0特殊混合型の有限七軌道被覆: `S61_R0_MIXED_ORBIT_REDUCTION_20261101`。依存と独立レビューの範囲は各ノートに記載。
- 四つの混合曲線型の普遍最高階数上界: `S61_R0_MIXED_CURVE_UNIFORM_RANK_20261101`。依存と独立レビューの範囲は各ノートに記載。
- 四つの混合曲線軌道の全一般β除去: `S61_R0_MIXED_CURVE_EXHAUSTION_20261101`。依存と独立レビューの範囲は各ノートに記載。
- 混合曲面T7の四重構造をもつ局所開: `S61_R0_MIXED_SURFACE_T7_FORMAL_20261101`。依存と独立レビューの範囲は各ノートに記載。
- M1二基線の全有限部分の除去: `S61_MINIMAL_M1_BASELINE_REMOVAL_20261101`。依存と独立レビューの範囲は各ノートに記載。
- M0二基線の全有限部分の除去: `S61_MINIMAL_M0_BASELINE_REMOVAL_20261101`。依存と独立レビューの範囲は各ノートに記載。
- 六つの半単純最小型尾の全有限部分の除去: `S61_MINIMAL_SEMISIMPLE_TAILS_REMOVAL_20261101`。依存と独立レビューの範囲は各ノートに記載。
- 二平面型の深い十七次元階数降下層の自由65開: `S61_BINARY_DEEPER_FORMAL_20261101`。依存と独立レビューの範囲は各ノートに記載。
- K0上のB階数七以下の全四素層: `S61_BINARY_DEEPER_REDUCTION_20261101`。依存と独立レビューの範囲は各ノートに記載。
- `S61_BINARY_MIXED_CONNECTION_20261101` は横断上界67までの途中結果であり、全局所環の結論ではない。

## v40: T7全軌道・最小有理尾全有限部分と境界局所構造

- 混合曲面T7軌道の全一般β除去: `S61_R0_MIXED_SURFACE_T7_EXHAUSTION_20261101`。証明・係数証明書・独立レビューの範囲を参照。
- 混合曲面T5の一般点付近を次元比較で除去: `S61_R0_MIXED_SURFACE_T5_OPEN_20261101`。証明・係数証明書・独立レビューの範囲を参照。
- 旗型の全b0b1非零開を除去し二不変因子へ還元: `S61_FLAG_UNIFORM_OPEN_20261101`。証明・係数証明書・独立レビューの範囲を参照。
- 二平面・混合尾接続点の随伴次数環と非零平方零t: `S61_BINARY_MIXED_FORMAL_20261101`。証明・係数証明書・独立レビューの範囲を参照。
- 二つの最小三重尾の全有限部分の除去: `S61_MINIMAL_TRIPLE_TAILS_REMOVAL_20261101`。証明・係数証明書・独立レビューの範囲を参照。
- 六つの最小純三次尾の全有限部分の除去: `S61_MINIMAL_PURE_CUBIC_TAILS_REMOVAL_20261101`。証明・係数証明書・独立レビューの範囲を参照。
- 四つの最小混合尾の全有限部分と射影端点の除去: `S61_MINIMAL_MIXED_TAILS_REMOVAL_20261101`。証明・係数証明書・独立レビューの範囲を参照。
- `S61_BINARY_MIXED_RELATIVE_DOUBLE_LINE_20261101` は全H9の相対二重直線を与える補助結果で、相対全環のt²0はこの補題単独では主張しない。
- 二平面・混合16D層の全一般β除去: `S61_BINARY_MIXED_STRATUM_REMOVAL_20261101`。相対tf特殊化とambient忠実平坦降下を別に証明。

## v41: 主曲線全有限部分・二平面境界・旗接続自由70

- M0主曲線c=±1のJordan軌道を全一般βで除去: `S61_MINIMAL_MAIN_JORDAN_REMOVAL_20261101`。レビュー範囲と未解決部分を各ノートに明記。
- 最小M0/M1主曲線の全有限部分を除去: `S61_MINIMAL_MAIN_REMOVAL_20261101`。レビュー範囲と未解決部分を各ノートに明記。
- K=0上のB階数六以下を一般βから除外: `S61_BINARY_DEEPER_RANK_SIX_EXCLUSION_20261101`。レビュー範囲と未解決部分を各ノートに明記。
- K=0二平面族の残余線形階数を二場合へ還元: `S61_BINARY_HYPERSURFACE_LOWER_RANK_REDUCTION_20261101`。レビュー範囲と未解決部分を各ノートに明記。
- 深い二平面chartのf=0全境界二層を除去: `S61_BINARY_DEEPER_BOUNDARY_REMOVAL_20261101`。レビュー範囲と未解決部分を各ノートに明記。
- 深い二平面の二つのpivot因子を除去: `S61_BINARY_DEEPER_PIVOT_DIVISORS_20261101`。レビュー範囲と未解決部分を各ノートに明記。
- 実際の有限体截断環から標数零形式長への上界移送: `ADIC_MODULAR_LENGTH_BOUND_20261101`。レビュー範囲と未解決部分を各ノートに明記。
- 旗・binary接続点の自由70と非零平方零t: `S61_FLAG_BINARY_CONNECTION_FORMAL_20261101`。レビュー範囲と未解決部分を各ノートに明記。

## v42: T2全層除去・深部有理chart・T5自由23

- 等階数特殊化での有限底torsionと飽和境界: `FINITE_BASE_TORSION_SPECIALIZATION_20261101`。適用範囲と未解決部分を各ノートに明記。
- T2接続指定点の随伴次数環と平方零t: `S61_BINARY_T2_CONNECTION_20261101`。適用範囲と未解決部分を各ノートに明記。
- T2接続16次元層を全一般βで除去: `S61_BINARY_T2_RELATIVE_DOUBLE_LINE_20261101`。適用範囲と未解決部分を各ノートに明記。
- 深部二平面族を七変数有理切断へ還元: `S61_BINARY_DEEPER_SEVEN_SLICE_20261101`。適用範囲と未解決部分を各ノートに明記。
- 深部の普遍三角消去と追加二座標例外の除去: `S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101`。適用範囲と未解決部分を各ノートに明記。
- T5普遍法線長23下界と指定点の自由23: `S61_R0_T5_FORMAL_RANK23_20261101`。適用範囲と未解決部分を各ノートに明記。

## v43: 良い還元・端点B除去・普遍二次rank24

- 二平面普遍法線下界65の標数32003での良い還元: `S61_BINARY_PRIME_SPREADING_20261101`。適用範囲と未解決部分を各ノートに明記。
- 半単純尾の二平面端点Bを全一般βで除去: `S61_MINIMAL_SEMI_B_ENDPOINT_20261101`。適用範囲と未解決部分を各ノートに明記。
- 深部二平面族の全一般β点で二次階数24: `S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101`。適用範囲と未解決部分を各ノートに明記。
- 旗b1=0指定点の全135式の正確六次jet: `S61_FLAG_BOUNDARY_EXACT_JETS_20261101`。適用範囲と未解決部分を各ノートに明記。

## v44: 端点A・純三次binary端点・T5普遍開・旗QQ六次

- 半単純尾の二平面端点Aを全一般βで除去: `S61_MINIMAL_SEMI_A_ENDPOINT_20261101`。適用範囲と未解決部分を各ノートに明記。
- 純三次尾の非平方二平面端点二表示を全一般βで除去: `S61_MINIMAL_PURE_BINARY_ENDPOINT_20261101`。適用範囲と未解決部分を各ノートに明記。
- T5の明示的開全点で自由階数23: `S61_R0_T5_OPEN_RANK23_20261101`。適用範囲と未解決部分を各ノートに明記。
- T5の二次pivot例外二層を全て除去: `S61_R0_T5_QUADRATIC_ARTIFACTS_20261101`。適用範囲と未解決部分を各ノートに明記。
- 旗境界のQQ全六次截断idealと初期部分の構造: `S61_FLAG_BOUNDARY_RATIONAL_INITIAL_20261101`。適用範囲と未解決部分を各ノートに明記。
- T5の全一般β軌道を飽和除去: `S61_R0_T5_ZERO_L_EXHAUSTION_20261101`。固定βの同じ曲面と相対smooth降下、四場合の全被覆を確認。

## v45: 三重尾・rank-one端点・旗固定点と有限正常化

- 三重尾の二平面端点を全一般βで除去: `S61_MINIMAL_TRIPLE_ENDPOINT_20261101`。適用範囲と未解決部分を各ノートに明記。
- rank-one A1族の全射影端点を一般βで除去: `S61_MINIMAL_RANKONE_ENDPOINT_20261101`。適用範囲と未解決部分を各ノートに明記。
- 旗境界の疎な固定β点で非零平方零t: `S61_FLAG_BOUNDARY_ORBIT_SPLIT_FORMAL_20261101`。適用範囲と未解決部分を各ノートに明記。
- 旗b1=0境界の有限正常化と線形例外の全域判定: `S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101`。適用範囲と未解決部分を各ノートに明記。

## v46: 旗境界の低次・端点表・深部三次

- 残る十二旗端点表示を四つの優越十六次元層へ整理: `S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101`。適用範囲と未解決部分を各ノートに明記。
- 深部binaryの二つの真三次制約が二core変数へ降りる: `S61_BINARY_DEEPER_CUBIC_NUISANCE_20261101`。適用範囲と未解決部分を各ノートに明記。
- 八変数旗chartの全真三次を二パラメータへ簡約: `S61_FLAG_BOUNDARY_EIGHT_CUBIC_20261101`。適用範囲と未解決部分を各ノートに明記。
- 旗八変数chartの二次rank9を線形枠開全体で被覆: `S61_FLAG_BOUNDARY_QUADRATIC_COVER_20261101`。適用範囲と未解決部分を各ノートに明記。
- b12=0旗七変数chartの二次・三次を全域で確定: `S61_FLAG_ZERO12_LOW_JETS_20261101`。適用範囲と未解決部分を各ノートに明記。
- b1=0十八次元旗境界の非空ambient除去開: `S61_FLAG_BOUNDARY_RELATIVE_REMOVED_OPEN_20261101`。適用範囲と未解決部分を各ノートに明記。
- r=−91/400旗層の真三次と異なるHilbert関数: `S61_FLAG_RANKONE_SOCLE_CUBIC_20261101`。適用範囲と未解決部分を各ノートに明記。
- 自由M0旗端点の実在平方零法線三方向を持つ平坦商: `S61_MINIMAL_FREE_FLAG_FLAT4_PROBE_20261101`。適用範囲と未解決部分を各ノートに明記。
- 乗法rank1 core全点に固定βの実在T7曲線を構成: `S61_FLAG_T7_BOUNDARY_20261101`。厚化・飽和は未主張。

## v47: 旗境界の真四次・T7平坦商・深部四次

- b12=0旗七変数chart全域の真四次とtの位数: `S61_FLAG_ZERO12_QUARTIC_20261101`。適用範囲と未解決部分を各ノートに明記。
- T7接続coreのb8開全域で真四次を確定: `S61_FLAG_RANKONE_SOCLE_QUARTIC_20261101`。適用範囲と未解決部分を各ノートに明記。
- T7曲面を旗端点まで延ばす実在平坦四重商: `S61_FLAG_T7_FLAT4_20261101`。適用範囲と未解決部分を各ノートに明記。
- 深部binaryの共通core開上で横断長65上界: `S61_BINARY_DEEPER_QUARTIC_NUISANCE_20261101`。適用範囲と未解決部分を各ノートに明記。
- 自由M0旗端点の有限体十次初期idealを再現: `S61_MINIMAL_FREE_FLAG_INITIAL_PROBE_20261101`。適用範囲と未解決部分を各ノートに明記。
- 四つの最小型旗端点層の独立監査: `S61_MINIMAL_FLAG_ENDPOINT_REVIEW_20261101`。適用範囲と未解決部分を各ノートに明記。
- b12=0旗chartの実在C3と指定点の有限rank5: `S61_FLAG_ZERO12_RANK_FIVE_20261101`。指定点のrankと全域のsupport幾何を区別。
- 自由M0旗端点の全一般βを覆う五変数正常形: `S61_MINIMAL_FREE_FLAG_NORMAL_FORM_20261101`。三つの悪い像も処理した有限正常形。

## v48: 深部core開の除去・旗横断49・T7全core被覆

- T7接続coreの全一般β点に真四次を拡張: `S61_FLAG_RANKONE_SOCLE_ZERO8_QUARTIC_20261101`。適用範囲と未解決部分を各ノートに明記。
- T7旗端点の平坦四重商を全ordinary coreへ拡張: `S61_FLAG_T7_FLAT4_COVER_20261101`。適用範囲と未解決部分を各ノートに明記。
- 深部binary共通core開の全一般β点を飽和除去: `S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101`。適用範囲と未解決部分を各ノートに明記。
- 自由M0旗端点の五変数普遍二次: `S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUADRATIC_20261101`。適用範囲と未解決部分を各ノートに明記。
- 自由M0旗端点の真三次を定数化: `S61_MINIMAL_FREE_FLAG_UNIVERSAL_CUBIC_20261101`。適用範囲と未解決部分を各ノートに明記。
- b12=0旗chart全域で実横断族は自由49: `S61_FLAG_ZERO12_TRANSVERSE49_20261101`。適用範囲と未解決部分を各ノートに明記。
- 自由M0五変数正常形と低次の独立監査: `S61_MINIMAL_FREE_FLAG_LOW_REVIEW_20261101`。適用範囲と未解決部分を各ノートに明記。
- 自由M0五変数の真正四次を確定: `S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUARTIC_20261101`。有限四次の定理であり独立四次監査は未実施。

## v49: 二次core例外の除去・旗真四次・最小型真五次

- 通常八変数旗の真正四次とt次数: `S61_FLAG_BOUNDARY_EIGHT_QUARTIC_20261101`。適用範囲と未解決部分を各ノートに明記。
- 通常旗四次の一変数小行列式被覆: `S61_FLAG_BOUNDARY_EIGHT_QUARTIC_COVER_20261101`。適用範囲と未解決部分を各ノートに明記。
- b12=0指定点のQQ六次初期idealと埋込111: `S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101`。適用範囲と未解決部分を各ノートに明記。
- 自由M0五変数の真正五次を確定: `S61_MINIMAL_FREE_FLAG_UNIVERSAL_QUINTIC_20261101`。適用範囲と未解決部分を各ノートに明記。
- 自由M0四次の独立低次商監査: `S61_MINIMAL_FREE_FLAG_QUARTIC_REVIEW_20261101`。適用範囲と未解決部分を各ノートに明記。
- 深部binaryの残る二次minor八因子層を除去: `S61_BINARY_DEEPER_CORE_ARTIFACTS_20261101`。適用範囲と未解決部分を各ノートに明記。

## v50: 旗四次の全基底被覆・七次閉包・最小型真六次

- 通常旗四次の全高次基底例外を解消: `S61_FLAG_BOUNDARY_GLOBAL_QUARTIC_COVER_20261101`。適用範囲と未解決部分を各ノートに明記。
- 通常旗h=0の四次を補完: `S61_FLAG_BOUNDARY_ZERO8_QUARTIC_COVER_20261101`。適用範囲と未解決部分を各ノートに明記。
- b12零指定点の真正七次閉包: `S61_FLAG_ZERO12_SEVENTH_INITIAL_20261101`。適用範囲と未解決部分を各ノートに明記。
- 全初期代数と埋込自由性だけでは冪零性は従わない: `FILTERED_INITIAL_NILPOTENCE_OBSTRUCTION_20261101`。適用範囲と未解決部分を各ノートに明記。
- 自由M0の真正六次に有限切断上界: `S61_MINIMAL_FREE_FLAG_SIXTH_UPPER_BOUND_20261101`。適用範囲と未解決部分を各ノートに明記。
- 自由M0の真正六次を普遍的に確定: `S61_MINIMAL_FREE_FLAG_UNIVERSAL_SEXTIC_20261101`。適用範囲と未解決部分を各ノートに明記。

## v51: 旗三線形例外・真rank降下の普遍三次と有限rank

- 旗の三線形基底例外を四次まで補完: `S61_FLAG_LINEAR_ARTIFACTS_QUARTIC_20261101`。適用範囲と未解決部分を各ノートに明記。
- 真の旗rank降下で全七変数二次を確定: `S61_FLAG_TRUE_CORE_QUADRATIC_20261101`。適用範囲と未解決部分を各ノートに明記。
- 真の旗coreの全七変数真正三次: `S61_FLAG_TRUE_CORE_CUBIC_20261101`。適用範囲と未解決部分を各ノートに明記。
- 真の旗core指定点の有限基底rank5: `S61_FLAG_TRUE_CORE_RANK_FIVE_20261101`。適用範囲と未解決部分を各ノートに明記。
- 自由M0の七次有限上界: `S61_MINIMAL_FREE_FLAG_SEVENTH_UPPER_BOUND_20261101`。適用範囲と未解決部分を各ノートに明記。
- 自由M0の八次までの累積上界: `S61_MINIMAL_FREE_FLAG_EIGHTH_UPPER_BOUND_20261101`。適用範囲と未解決部分を各ノートに明記。
- 真の旗coreの全七変数真正四次: `S61_FLAG_TRUE_CORE_QUARTIC_20261101`。適用範囲と未解決部分を各ノートに明記。

## v52: AB全域正規性と実在の旗平方零商

- AB正規化の四次を完全確定: `AB_ORIGIN_DEGREE_FOUR_20261101`。適用範囲と未解決部分を各ノートに明記。
- ABの原点を含む全域正規性: `AB_ORIGIN_NORMALITY_20261102`。適用範囲と未解決部分を各ノートに明記。
- 全b0非零旗境界の実在閉平方零厚化: `S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101`。適用範囲と未解決部分を各ノートに明記。
- 深部binary七二次idealの積正規形: `S61_BINARY_DEEPER_PRODUCT_PAIR_20261102`。適用範囲と未解決部分を各ノートに明記。

## v53: AB有理特異点・T01零fiberの一般被約性・有限模型

- ABの有理特異点と全Hilbert級数: `AB_RATIONAL_SINGULARITIES_HILBERT_20261102`。適用範囲と未解決部分を各ノートに明記。
- 実際のT01零fiberの一般被約性: `TORIC_ZERO_FIBER_GENERIC_REDUCED_20261102`。適用範囲と未解決部分を各ノートに明記。
- zero12の全基底次数に対する有限打切り: `S61_FLAG_ZERO12_BASE_ADIC_BOUNDS_20261102`。適用範囲と未解決部分を各ノートに明記。
- 深部binaryの厳密有限weighted模型: `S61_BINARY_DEEPER_WEIGHTED_FINITE_MODEL_20261102`。適用範囲と未解決部分を各ノートに明記。

## v54: 深部二平面型・旗境界の除去と残る三core

- 深部二平面型の三次階数と十七曲線: `S61_BINARY_WEIGHTED_CUBIC_COVER_20261102`。適用範囲と未解決部分を各ノートに明記。
- 深部二平面型での四次階数欠損の有限性: `S61_BINARY_WEIGHTED_HORIZONTAL_SUPPORT_20261102`。適用範囲と未解決部分を各ノートに明記。
- 深部二平面型全層の一般β上の除去: `S61_BINARY_DEEPER_FULL_CORE_REMOVAL_20261102`。適用範囲と未解決部分を各ノートに明記。
- 相対冪零性の正重み収縮: `S61_BINARY_NUISANCE_CONTRACTION_20261102`。適用範囲と未解決部分を各ノートに明記。
- 隣接超平面層の八法線有限表示: `S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102`。適用範囲と未解決部分を各ノートに明記。
- 三core超平面層での全三次階数欠損: `S61_BINARY_HYPERPLANE_CUBIC_CODIMENSION_20261102`。適用範囲と未解決部分を各ノートに明記。
- 二変数関数体上の全dにわたる自由65: `S61_BINARY_HYPERPLANE_BFREE_HORIZONTAL_20261102`。適用範囲と未解決部分を各ノートに明記。
- 正規化H8の一般β上の残点の有限性: `S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_20261102`。適用範囲と未解決部分を各ノートに明記。
- 指定zero12の第一底商: `S61_FLAG_ZERO12_FIRST_BASE_QUOTIENT_20261102`。適用範囲と未解決部分を各ノートに明記。
- 実際の二変数有限核と低次作用: `S61_FLAG_ZERO12_FIRST_LINEAR_T_ACTION_20261102`。適用範囲と未解決部分を各ノートに明記。
- zero12旗境界全点の一般β除去: `S61_FLAG_ZERO12_RELATIVE_WEIGHT_REMOVAL_20261102`。適用範囲と未解決部分を各ノートに明記。
- 三値を除く通常旗境界全域の一般β除去: `S61_FLAG_ORDINARY_CORE_SUPPORT_20261102`。適用範囲と未解決部分を各ノートに明記。
- 残る三coreの固定β環: `S61_FLAG_REMAINING_CORE_WEIGHTED_POINTS_20261102`。適用範囲と未解決部分を各ノートに明記。
- 自由M0指定有限体点の底三次までの実商: `S61_MINIMAL_FREE_FLAG_THIRD_BASE_QUOTIENT_20261102`。適用範囲と未解決部分を各ノートに明記。
- ABの明示二次基底と三次関係の尽くし: `AB_EXPLICIT_QUADRICS_CUBIC_GENERATION_20261102`。適用範囲と未解決部分を各ノートに明記。

## v55: ABの二次生成・一般P25五点断面と旗の局所構造

- ABの素イデアルの明示二次生成: `AB_QUADRICS_QUARTIC_MULTIPLICATION_20261102`。適用範囲と未解決部分を各ノートに明記。
- 平面外の全射影解の法線方向からの復元: `S61_AFFINE_NORMAL_PROJECTION_20261102`。適用範囲と未解決部分を各ノートに明記。
- 一般P25切断でのスキーム論的五点の尽くし: `S61_AFFINE_TWENTYFIVE_PLANE_20261102`。適用範囲と未解決部分を各ノートに明記。
- 全非平方二平面型の階数34以下での残余有限性: `S61_BINARY_H11_RANK_DROP_RESIDUAL_FINITE_20261102`。適用範囲と未解決部分を各ノートに明記。
- 指定T7横断環の一次分解: `S61_FLAG_T7_WEIGHTED_PRIMARY_20261102`。適用範囲と未解決部分を各ノートに明記。

## v56: 実T01零ファイバー、分岐、一般アフィン支持

- 実T01零ファイバーの全域での被約性: `TORIC_ZERO_FIBER_COFRAME_QUADRICS_20261102`。適用範囲と未解決部分を各ノートに明記。
- 判別因子の一般点での実T01のファイバー: `TORIC_GENERIC_DISCRIMINANT_AFFINE_FIBER_20261102`。適用範囲と未解決部分を各ノートに明記。
- 法線20変数上の有限表示と射影26平面: `S61_AFFINE_FINITE_NORMAL_MODULE_20261102`。適用範囲と未解決部分を各ノートに明記。
- 未知アフィン部分を消す具体的な行列式: `S61_AFFINE_RESIDUAL_ANNIHILATOR_20261102`。適用範囲と未解決部分を各ノートに明記。
- T7二切断の実相対torsion: `S61_FLAG_T7_TRANSVERSE_BASE_TORSION_20261102`。適用範囲と未解決部分を各ノートに明記。

## v57: 一般アフィン次元の上界・T7境界と判別因子の五平面

- 全型での真の無限遠境界の次元評価: `S61_AFFINE_BOUNDARY_DIMENSION_20261102`。適用範囲と未解決部分を各ノートに明記。
- 全一般アフィンファイバーの次元上界2: `S61_AFFINE_SURFACE_BOUND_20261102`。適用範囲と未解決部分を各ノートに明記。
- 元の二座標u6,u42への有限射: `S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102`。適用範囲と未解決部分を各ノートに明記。
- 射影零ファイバーと次数18の極方向写像: `TORIC_PROJECTIVE_ZERO_FIBER_POLAR_IMAGE_20261102`。適用範囲と未解決部分を各ノートに明記。
- 判別因子の線形五平面と極方向の逆算: `TORIC_DISCRIMINANT_LINEAR_RULING_20261102`。適用範囲と未解決部分を各ノートに明記。
- T7旗core全層の一般β上の飽和除去: `S61_T7_CRITICAL_AMBIENT_TRANSVERSE_20261102`。適用範囲と未解決部分を各ノートに明記。

## v58: 通常旗の全域除去と一般アフィン曲面の次数制約

- T5通常旗の全域除去: `S61_FLAG_T5_CORE_REMOVAL_20261102`。適用範囲と未解決部分を各ノートに明記。
- 正重み底を残すr=2の全域除去: `S61_FLAG_R2_POSITIVE_BASE_REMOVAL_20261102`。適用範囲と未解決部分を各ノートに明記。
- 指定因子が重なる平方旗の残余有限性: `S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102`。適用範囲と未解決部分を各ノートに明記。
- 関数体上の境界の剰余次数: `S61_GENERIC_BOUNDARY_INDEX_20261102`。適用範囲と未解決部分を各ノートに明記。
- 仮想アフィン曲面の次数下界15: `S61_AFFINE_SURFACE_DEGREE_OBSTRUCTION_20261102`。適用範囲と未解決部分を各ノートに明記。
- 固定二座標上での166生成表示: `S61_AFFINE_FIXED_PAIR_LOCALIZATION_20261102`。適用範囲と未解決部分を各ノートに明記。
- 五平面と節点復元の独立点検: `TORIC_DISCRIMINANT_LINEAR_RULING_REVIEW_20261102`。適用範囲と未解決部分を各ノートに明記。

## v59: 非平方零対の除去と有限境界点への帰着

- c=f=0の非平方二平面型全層の除去: `S61_BINARY_ZERO_PAIR_STRATUM_REMOVAL_20261102`。適用範囲と未解決部分を各ノートに明記。
- 全境界曲線を有限個の指定点へ帰着: `S61_AFFINE_BOUNDARY_FINITE_ANCHORS_20261102`。適用範囲と未解決部分を各ノートに明記。
- 有限集合の局所商による曲面次数の上界: `S61_AFFINE_ANCHOR_LENGTH_BOUND_20261102`。適用範囲と未解決部分を各ノートに明記。
- 指定平方旗のb1=0交叉の全除去: `S61_MARKED_SIGMA_INTERSECTION_REMOVAL_20261102`。適用範囲と未解決部分を各ノートに明記。
- 指定平方旗の射影coreと実アフィン直線: `S61_MARKED_SIGMA_PROJECTIVE_CORE_20261102`。適用範囲と未解決部分を各ノートに明記。
- 166生成表示の一変数表示: `S61_AFFINE_FIXED_PAIR_UNIVARIATE_20261102`。適用範囲と未解決部分を各ノートに明記。

## v60: 有限射影候補、境界の剰余次数と平方旗因子

- 非平方零対の独立再現: `S61_BINARY_ZERO_PAIR_STRATUM_REMOVAL_REVIEW_20261102`。適用範囲と未解決部分を各ノートに明記。
- 真の境界が会う有限点の剰余次数: `S61_AFFINE_ANCHOR_RESIDUE_DEGREE_20261102`。適用範囲と未解決部分を各ノートに明記。
- 指定平方旗の純三次端点因子の除去: `S61_MARKED_SIGMA_PURE_ENDPOINT_REMOVAL_20261102`。適用範囲と未解決部分を各ノートに明記。
- 指定平方旗の剰余次数15候補の全除去: `S61_MARKED_SIGMA_ZERO_H_REMOVAL_20261102`。適用範囲と未解決部分を各ノートに明記。
- 指定平方旗の試験直線全域の方程式: `S61_MARKED_SIGMA_PENCIL_MODEL_20261102`。適用範囲と未解決部分を各ノートに明記。
- 元135式の三次元基本サイクル: `S61_PROJECTIVE_TOP_CYCLE_20261102`。適用範囲と未解決部分を各ノートに明記。
- 座標二つを固定した有限な射影候補: `S61_FIXED_PAIR_SINGLE_COLON_20261102`。適用範囲と未解決部分を各ノートに明記。
- 全20法線方向での斉次化変数の有限作用: `S61_AFFINE_HOMOGENIZER_FINITE_ACTION_20261102`。適用範囲と未解決部分を各ノートに明記。
