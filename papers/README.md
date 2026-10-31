# 公表用原稿の構成案と証明資料

**状態：v28-s61、2026-10-31。一次部分の研究草稿を統合更新。**

[一次ジェット還元の日本語原稿](s61_first_jet_reduction_ja.tex)（[PDF](s61_first_jet_reduction_ja.pdf)、22頁）は、
primitive gauge、T01/T39/T40の優越性と全S61のétale点、T01の一般次数3、
少なくとも五つの一般延長、トーリック代数模型と無限遠の旗、階数1固定fibre、
CA/CBの共通像を扱う。MPS40の全図式と内部IDも一次資料へ照合済み。
内部研究草稿であり、新規性の網羅的調査と投稿準備は未完。
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
