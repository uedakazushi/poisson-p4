# 研究状態 — 全分類の共通復元と実持上げ

2026-11-09 / v68-marked-reconstruction

active queueは `research/conceptual/GLOBAL_PLAN.md`。旧classification/PLAN.jsonは構築完了の記録である。
現行分類は147成分（rank4:141、rank2:6）、1158採録主張を主張する。
その元の全支持証明は `research/classification/GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF.md`、
再現入口は `research/classification/REPRODUCE_EXACT147.md`。今回全147を独立検証したとはしない。

単純Hと重複Hの入口の重複を、共有された標識超平面・modular原始関数・共鳴の証明へ置換した。
第三の二H+三次模型、定数modular全源の正逆と加法定数も保持する。
TRIPLE→Ga9では次数制約から補正を導き、指定因子族の二次持上げに補正が必要なことを示した。
別の可換場の分裂族では全高次補正が零である。任意方向の無障害性を主張しない。
主稿の全既存証明を保持し、成分台帳・主張台帳・旧証明書は変更していない。

対応表はactual因子記録とMPSの一般源から作成し、A113/B10/C18/R2六。
両原稿は共通TeX二本を読む。新規性未確認、全source帰属と成分性の旧依存は保持する。
検算は `code/check_marked_reconstruction.py` と既存 `code/check_global_structure.py`。
TeXは `scripts/build_structural_manuscripts.sh` で両方二回組版する。
実結果は `reports/marked_reconstruction_validation.json`。

次は超平面なしQ+Cのcanonical symmetryと全31代表を閉形式復元へ統合する。
続いてCの導手/正規化とrank2支持、極大閉像を接続する。H+2Q局所補題は必須前提ではない。
GitHub保存、bundle生成/検証、Library保存は別状態として実際の結果だけを報告する。
