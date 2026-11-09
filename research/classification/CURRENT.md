# 現在地 — 完全分類完成報告と公表前の独立査読

2026-11-09 / v73-prepublication-review

FINAL_TARGETのM1–M4、C1–C2、P1、R1、S1を満たしたと原稿と保存記録は報告する。主稿の独立した公表前査読は別の検討であり、現段階で完了したとは報告しない。
最終仕様は `research/conceptual/FINAL_TARGET.md`。通常終了条件は維持する。
全147の証明を今回独立に再証明したという意味ではなく、下記の既存証明・一次資料・正確な証明書を指定の射程で合成する。

| 条件 | 現行の達成証拠 |
|---|---|
| M1 全点被覆 | 主稿prop:global-exact147。A/B/Cの全source・正逆・零/分母補集合とprop:irreducible-nonnormal-all-pointwiseを保持。rank四閉包を先に覆い、rank_two_support_reconstruction.texのprop:rank2-supportへ戻す。安定な五切断、非安定・因子付きの全枝、toric/Aff32の実変形を統合 |
| M2 極大性・相違 | maximal_image_synthesis.tex。全周囲の挟み撃ち、真正形式切片、有限閉和からの極大性を合成。全既存142との必要比較、最後五像の相違、標識/Galois/frame同定、特殊底を処理。cor:reduced-divisor-component-imagesの像一致禁止を包含判定と混同しない |
| M3 具体的記述 | component_family_catalogue.texに全147 IDを一対一で定義。構成・成分性・相違を別参照にし、conic+cubic六型、二H+cubicの追加三型、MU5四型、nil4 quartic、NEW25等の実命題を照合。最後五族は全括弧を明記。51族の有限有理演算とS61復元、NEW25全多項式源、各全sourceの枝・開条件・次元・階数・因子/foliation・証明を指定。主稿のcatalogue:IDと案内CSVを一致させた |
| M4 依存閉鎖 | 主稿・必要別稿の証明と既存有限証明書を再利用。新たにKaledin Theorem1.5/Cor1.6、CM本文Theorem4.1/Prop3.9、LPT Theorem6.2、MPS Theorem1.4/6.16/Cor6.23/Example6.24の仮定と適用を一次資料で確認。既存の正規因子認識の全証明を同時配布 |
| C1/C2 概念的統合 | marked H/modular原始関数→quadric cofactor/磁気商→余法線/全normal項→正規化/全B/A降下/導手pole/全scroll留数→実ambient変形→rank二支持→極大閉像。共通定理がdegree5/6と二つの実T40曲線等の重複を置換。有限認識・比較・残るJacobiを保持 |
| P1 公表原稿 | 主稿418頁、共通構造稿41頁、必要別稿の正規実因子認識91頁・一次還元118頁。全て二回組版し、box/未解決/重複参照ゼロ。改訂箇所を画像確認。旧global未完記述を訂正し、主証明から独立なS61方法論を付録に保持。優先権は未確認 |
| R1 再現一式 | required_proof_dependency_closure.json：223 primary claim、記録済み係数/コード入力の不一致・欠落ゼロ。96凍結入力を正確回収し、必要legacyを個別ファイル化。最新五分割ZIPから3579ファイルと相対参照2個を復元し、全SHA256と論理パスを照合。final_bundle_validation.jsonに記録 |
| S1 保存 | 現行の原稿・必要入力はsource commit a5f82342870b12583da0d0ca40bb5ba96c742b0eをmainで照合。Libraryに四PDFと五volumeを保存し、全access/hash及びPage sequence4の九リンクを読み戻して確認。final_library_persistence.jsonに記録 |

各検証の射程は `reports/final_manuscript_validation.json`、`reports/component_catalogue_integration_validation.json`、
`reports/rank_two_integration_validation.json`、`reports/normalization_residue_validation.json`、
`reports/dependency_reconstruction_checks.json` と元の各証明書にある。
新37正規化・28rank二検算、51演算列の正確標本、六rank二成分の元行列の再実行、
隔離したmixed rank3の全350/正確行列/二次障害−140の再実行を行った。
標本・integrity・組版を一般証明や全147の独立自動認証とは呼ばない。
font/scaling/hyperrefの既存警告は残る。

## 今回の追加フェーズ：公表前の独立査読

**主分類の構築完了報告をそのまま査読済みとは扱わない。**
研究方針は RESEARCH_POLICY.md、査読計画は research/publication/REFEREE_PLAN.md。
中心課題は thm:reduced-fibre-linear とその閉因子像の系、安定rank2のRat113/122全源、
最後の五成分と既存142との全非包含、51族の一般パラメータS61表示、
223 primary claim・96凍結入力の重要依存、五分割現行bundleの隔離復元である。
今回はこれらの独立査読を実行したとの主張ではなく、最も影響が大きい検査順を確定した。
既存の数学的誤りを今回発見したという主張もしない。
原稿と査読の範囲を分け、確認で欠陥が見つかれば関係する主張全体を訂正する。

## 保存確認

主稿は418頁、共通構造稿41頁、正規実因子認識91頁、一次還元118頁。
本文の構成・成分性・相違の参照を別々に照合し、最後五族の全括弧式も本文へ明記した。
reports/catalogue_reference_review.jsonに対応修正と五つの全パラメータJacobi/Pf検算を固定する。
族・次元・階数・147という結論は変更していない。

現行source commitはa5f82342870b12583da0d0ca40bb5ba96c742b0e。
最終bundleのbyte hash、復元、Libraryの参照と読み戻しは
reports/final_bundle_validation.json及びreports/final_library_persistence.jsonにある。
資料ページ：https://chatgpt.com/space/page_98c9971dc7c48191a2c4b7624fd99fa0

Libraryには全九ファイルの取得可能性を確認し、既存ファイルは削除していない。
自己参照hash循環を避け、bundleのsource commitと本保存receiptのcommitは分ける。
投稿又は新たな外部公開は行っていない。
正確な全PGL軌道、全交叉・全nil厚み、及び別稿の精密局所代数の残問題は主分類の終了条件外である。
