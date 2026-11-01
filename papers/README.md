# 公表用原稿の構成案と証明資料

**状態：v45-s61、2026-11-01。零一次部分の被約台AB、ABのincidenceと特異点集合、最小型無限遠の完成局所環、二平面型の横断Artin代数を統合。**

[一次ジェット還元の日本語原稿](s61_first_jet_reduction_ja.tex)（[PDF・41ページ](s61_first_jet_reduction_ja.pdf)）は、
primitive gauge、T01/T39/T40の優越性と全S61のétale点、楕円型五変数模型、T01の一般次数3、
少なくとも五つの一般延長、トーリック代数模型と無限遠の旗、階数1固定fibre、
CA/CBの共通像、全模型の普遍的な重み6帰結を扱う。MPS40の全図式と内部IDも一次資料へ照合済み。
v31では標準曲線の既約被約判別式、Jacobianの非零二捩れとの同型、実際のT01のΔ非零上の有限エタール性、全域の飽和イデアルと欠落境界を追加した。
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
