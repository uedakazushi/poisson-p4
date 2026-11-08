# 公表前監査の現在地と再開手順

2026-11-09 / v1。対象commit: ffc178bff0a69fb0edcc99a4d8d7de98c8a2a222。

## 結論の確実性

現行主稿は147成分の完全分類を主張する。今回その数を変更する理由や反例は得ていない。
ただし全主張を独立に検証したわけではない。次の表は今回実際に確認した範囲であり、過去のPASSの転記ではない。

| 対象 | 今回の状態 | 証明・証拠 |
|---|---|---|
| 有限閉既約被覆と相互非包含から成分分類 | 証明を確認・明文化 | verification_note_ja.tex, lem:cover |
| rank4全点被覆の後でrank2閉包項を戻す論理 | 証明を確認・明文化 | 同 lem:rank |
| TRIPLE→Ga9の表示された一族 | 独立に再実装し81項目を正確計算 | code/verify_final_bridge.py, reports/final_bridge.json |
| 非線形補助写像と真の射影線形帰属の区別 | 両方の写像・逆・全括弧を確認 | 同 prop:family, prop:ga9, prop:closure |
| 低次数不変曲面の特殊化機構 | 仮定付き一般補題を証明 | 同 lem:specialize |
| Corrêa–Muniz v5 Theorem4.1、Corollary D | 一次HTML本文の主張と適用範囲を確認 | 文献欄 |
| MPS v2 Theorem1.4、Example6.24 | 一次HTML本文の主張と適用範囲を確認 | 文献欄 |
| 147全成分の全点被覆・個別成分性・相違 | 今回の独立検証は未完 | 下の重点課題 |

81はコード内の恒等式・次数・非消失等の検査項目数であり、81個の独立な数学定理を意味しない。
旧150/350/126形式を再生成したという報告ではない。独立なh-Casimir斉次二次表示でJacobiと退化因子を検証した。
根拠と証明はLuaLaTeXノートに置き、ここでは重複しない。

## 優先する原証明と着手方法

**P1：全rank4被覆。** `research/classification/KNOWN142_SEVEN_CANDIDATE_WHOLE_POISSON_COVER_AUDIT.md` の各行を終端定理へ辿る。
特に `papers/normal_actual_divisor_recognition_ja.tex` と `research/classification/IRREDUCIBLE_NONNORMAL_QUINTIC_WHOLE_POINTWISE_COVERAGE.md` を優先する。
各行に「原仮定／全パラメータ源／逆写像／零・単元補集合／帰属先／主稿の命題」を対応させる。
上位の合成監査を読むだけでは終了しない。旧sourceや証明書のハッシュを変更せず、誤りがあれば依存する現行主張を訂正する。

**P2：rank2支持。** `research/classification/RANK_TWO_SUPPORT.md` と `RANK2_STABLE_REVIEW.md` から安定分岐の原証明へ。
次数2のunmixed曲面の分類、reflexive extension、五切断が全切断を張らない場合、発散条件を重点的に追う。
今回確認したCorrêa–Muniz Theorem4.1は安定分岐の代数的可積分性とuniruled葉を述べるが、それだけで二つの固定Jacobian族への所属は従わない。
因子付き分岐、Aff32の全源への稠密性、rank4への変形は別の証明義務である。

**P3：最後の五成分。** `T37_SIX_TARGET_LEAF_OBSTRUCTIONS_SEPARATE_REVIEW.md`（research/classification/内）のT37稠密軌道・全residue恒等式・対象側葉次数を原式で確認する。
今回の特殊化補題はその一部だけである。実D自身の滑らかな単純因子点という仮定を保持する。
`FIVE_CANDIDATES_FULL_COMPONENTS_AND_ONE_ORBIT_FRONTIER.md` の全既知親への非包含と候補相互非包含も別に確認する。

**P4：最後の閉包帰属の残入力。** 表示されたTRIPLE標準形のGa9への退化は今回確認済み。
`CONSTANT_MODULAR_TRIPLE_GA9_FULL_SOURCE_SEPARATE_REVIEW.md` の元全四源から一軌道への帰着、および `GA9_CONSTANT_COMPONENT.md` の全十源の成分性はまだ独立再検証していない。
一標本の成功を全sourceの成功へ昇格させない。

## 出版上の整理

主稿の証明を削除・要約しない。まず全依存を追跡し、その後、主定理・全点被覆・個別成分・計算証明の順に再構成する。
368頁の分割は未実施。現在の分類主稿・PDF・COMPONENTS.json・CLAIMS.jsonを今回変更していない。
既存papers/README.mdやPROOF_ROUTE.mdには歴史的な記述が残るので、現在の優先順位の入口としない。
旧証明資産の整理は依存確認後に行い、先に削除して軽量化しない。

## 一次文献

- Maurício Corrêa and Alan Muniz, Holomorphic foliations of degree two and arbitrary dimension, arXiv:2207.12880v5, 2026-01-19. Theorem 4.1、Corollary D。https://arxiv.org/html/2207.12880v5
- Mykola Matviichuk, Brent Pym and Travis Schedler, A local Torelli theorem for log symplectic manifolds, arXiv:2010.08692v2, 2021-07-15. Theorem 1.4、Example 6.24。https://arxiv.org/html/2010.08692v2

CMのrank2部分空間の成分数と、全Poisson空間のgeneric-rank2成分数は違う。MPSの40はSNC locusと交わる成分の分類である。
上記は2026-11-09の本文確認範囲。文献全体の証明再検証、新規性の徹底調査完了、外部査読済みという意味ではない。

## 再現と配布の制約

検算はREADMEの二コマンド。検証ノートの組版はLuaLaTeX二回、補助ファイルはwork/へ出力する。
今回の環境ではGitHubのテキスト読取とGit data APIによる更新は可能だが、主稿・PDFのコンテナへの直接ダウンロードに失敗した。
そのためLibrary用の今回の保存物は自己完結した最終接続監査パッケージであり、147成分分類全体の研究継続bundleを置き換えない。
完全な現行bundleの作成は、主稿・必要証明・証明書の依存閉包を取得できる環境で行う未完作業として残す。
