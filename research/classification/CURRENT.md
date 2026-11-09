# 現在地 — 共通復元と実持上げを主稿へ統合

2026-11-09 / v68-marked-reconstruction

既存分類稿・台帳は147成分（rank4:141、rank2:6）を主張する。今回その全証明を独立再確認したとはしない。
主分類稿の全既存証明を保持し、共通TeXを両原稿に読み込ませた。

- marked_hyperplane_reconstruction.tex：単純H+四次・重複Hの標識、modular、共鳴と最高部次数分岐を共通化。原始関数恒等式はPfaffianで割らず、定数modular全源と第三の二H+三次模型の正逆も保持。
- poisson_lifting_examples.tex：TRIPLE→Ga9の次数補正と二次持上げ、可換場による二重H分裂を通常の証明で説明。
- component_entrances.csv：全147 IDのactual一般元と出典。A113/B10/C18/R2六。閉包を分割する表ではない。

新しい恒等式・全表示族・正逆をSymPyで検算し、旧29項目も再現した。
個別源の成分性・相違・全支持は旧証明に依存し、成分・主張台帳と元証明書は未変更。新規性は未確認。
次はGLOBAL_PLAN P2のQ+C、導手、rank2と極大閉像。H+2Qの局所補題は必須前提ではない。
検証結果と保存範囲は reports/marked_reconstruction_validation.json およびbundle manifest。
