**Current persistence rule (2026-11-10):** GitHub main only; Library emergency only after failed push. No GitHub Actions. The operative reproduction contract is REPRODUCE_GITHUB.md, and the ongoing incidence study is governed by research/subvariety_incidence/RESEARCH_POLICY.md.

# 完全分類の数学的構造と現在の研究課題

2026-11-09 / v74-publication-candidate
本書は分類の現行構造を短く説明する。現在の実行計画は research/publication/REFEREE_PLAN.md。
指定した査読は公表候補版として可と判定した。根拠と配布状態は research/publication/REFEREE_REPORT.md と classification/CURRENT.md。
最終分類の数学的な完成仕様は FINAL_TARGET.md。読解・再開はRESEARCH_POLICY.mdを最初に読む。

全点認識 → 標識付き復元 → Poisson変形 → 相異なる極大閉像、という原理で全147成分の完成を記録している。
一般rank4には三つの入口 A: actual Dに超平面因子がある（5Hを含む）、
B: 超平面因子がなく既約二次＋既約三次、C: 被約既約五次の正規／非正規がある。
一般rank2は飽和foliationと零因子から扱い、rank4被覆の閉包を先に決めてから支持定理へ戻す。

現行の共通TeXは、marked_hyperplane_reconstruction、
quadric_cubic_reconstruction、conductor_conormal_reconstruction、
normalization_residue_reconstruction、rank_two_support_reconstruction、
poisson_lifting_examples、maximal_image_synthesis、
component_family_catalogue（いずれもpapers/）。
主稿 papers/poisson_p4_classification_ja.tex はこれらを利用する。
非正規五次の正規化・導手・全sourceと、最後五成分の非包含の原証明は
research/classification/GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF.md
および必要な元証明へ接続する。

全bundleの隔離検証とmain/Library保存が完了し、現在の目標を終了した。
新しい証拠又は具体的な誤りがある場合だけ、関連する証明依存を再開する。
確認済みでない箇所を「既に監査済み」に変更しない。
数学的な欠陥を見つけた場合だけ、新しい補題・一般化・訂正を開始する。
