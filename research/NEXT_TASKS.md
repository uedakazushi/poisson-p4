# 現行研究課題 — 非Cohen--MacaulayなPoisson成分交叉の局所構造

2026-11-10 / v11。GitHub mainのみ正本、TeX・計算コード・証明書を追跡。生成PDF非追跡、Actions不使用、Libraryは実際のpush失敗時のみ緊急退避。

**現行定理**：完全等方Log点p=[P0]の周囲PoissonスキームのZariski接空間30次元・局所次元28（分類の完全性を仮定しない）。同点はLog1112の28次元rank2成分とR4_MODULAR_RANK1_LONG_NONCONSTANT_25の25次元rank4成分との交点である。後者はPoisson曲線P0+aK1の実PGL5源同定により確認。二つの既約成分の次元が異なるため、局所環は非等次元・非CM・非Gorenstein・非局所完全交叉。結論の成分性入力は現行分類の証明に依存する。ambient非被約性は未証明。

**最優先**：他のrank4成分もpを通るか、二つの既知成分のscheme論的交叉と全接錐を定めること。既知の二つの障害式9*t21*t29、-4*t25*t27に続く独立二次式を計算し、分類カタログの他成分の閉包源と比較する。次に固定Qファイバーのac障害、ambient局所環の被約性を研究。論文はpapers/quadric_log_component_intersection.texとpapers/rank_four_long_component_membership.tex。計算はcode/subvariety_incidenceの新しい11番目のproducer。

再開：RESEARCH_POLICY.md → research/subvariety_incidence/RESEARCH_POLICY.md → NOTES.md → 現行TeX。手動再現：bash scripts/verify_github_source.sh incidence。
