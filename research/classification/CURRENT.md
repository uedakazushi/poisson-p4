# 現在地 — 147成分分類の公表前監査

2026-11-09 / v65-publication-review

分類構築の結論は正確147成分（一般rank4が141、rank2が6）。PGL5商前の射影Jacobiスキームの被約既約成分を扱う。
現行主稿は `papers/poisson_p4_classification_ja.tex` と同名PDF、既存ビルド報告は368頁・1158採録主張。
最終分類証明は `GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF.md`、再現入口は `REPRODUCE_EXACT147.md`。
これらの原稿・台帳と旧証明書を今回変更していない。

ユーザーの指示で現在の主作業を公表前監査へ移した。入口は `../publication/REVIEW.md`。
分類構築のPLAN.jsonのdoneと、今回の独立監査・外部査読の完了は異なる。旧READMEの124成分下界を現在地としない。

今回実行した独立検算はTRIPLE→Ga9明示族の81項目。積括弧から再構成し、全逆写像、全六括弧、斉次Jacobi、実Pfaffian、Ga9 potential、e=3,12を確認した。
有限閉被覆の極大成分判定、rank4閉包を先に扱う接続、不変曲面の特殊化は検証ノートで証明した。
全TRIPLE源の一軌道への帰着、Ga9成分性、正規・非正規五次、安定rank2、T37の源固有の計算は今回の検算範囲外。
主分類の具体的な反例を発見したという記録ではないが、全147成分の独立検証完了とも扱わない。

優先課題は原証明の全支持被覆と安定rank2/T37の接続の点検、その後に主稿の説明を再構成すること。
全軌道分類・全非被約厚み・精密nil/H0/T2等は引き続き対象外。
証明の詳細、仮定、確認範囲と次の着手方法は `../publication/REVIEW.md` を参照。
