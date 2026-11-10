## 現行v11研究入口（2026-11-10）

Poisson部分多様体incidenceの最新版は`research/subvariety_incidence/RESEARCH_POLICY.md` → `NOTES.md` → `papers/poisson_subvariety_incidence_ja.tex`。完全等方Log点p=[P0]の局所次元は28、接空間は30次元。pは28次元rank2 Log1112成分と25次元rank4長鎖非定数g成分（R4_MODULAR_RANK1_LONG_NONCONSTANT_25）の交点である。後者はPGL5の明示源同定による。既存の両成分性定理から局所環は非等次元・非CM。ambient非被約性および他のrank4成分は未決定。再現は`bash scripts/verify_github_source.sh incidence`（11 producer / 11 JSON）。GitHub main正本、PDF非追跡、Actions禁止、Libraryはpush失敗時のみ。

## GitHub正本規則
GitHub mainが研究の唯一の正本。再生成入口 REPRODUCE_GITHUB.md。新しい部分多様体incidenceは research/subvariety_incidence/RESEARCH_POLICY.md → NOTES.md → papers/poisson_subvariety_incidence_ja.tex。LibraryはGitHub push失敗時の緊急避難のみ、GitHub Actionsは禁止。下記は既存の147成分分類の完成状態。

# 研究状態の入口

2026-11-09 / v74-publication-candidate

現行の完全分類は147成分を主張し、最終成果物の完成・GitHubとLibraryの保存を報告している。
証明の全文を第三者的に独立査読したという意味ではない。
公表前の反証志向の査読は、指定した高リスク義務と有限依存の射程で完了した。
今回の定理読解・全51族の普遍Jacobi・原稿訂正は research/publication/REFEREE_REPORT.md。
65構文パスの読解と二runnerの欠落修復を反映し、公表候補版として可と判定した。
全bundleの隔離検証・main/Library保存も完了した。
唯一の実行状態は classification/CURRENT.md、保存receiptは reports/publication_candidate_delivery.json にある。

再開時：RESEARCH_POLICY.md → research/publication/REFEREE_PLAN.md →
research/classification/CURRENT.md → 検査対象の原証明。
構造理論の原稿・カタログ・全数値的・文献的入力は保持し、指示だけで数学の完成状態を捏造しない。
