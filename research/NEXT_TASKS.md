# 最優先 — Poisson部分多様体から射影四空間の成分分類を構造的に説明する

2026-11-11 / v18。再開：RESEARCH_POLICY.md → research/subvariety_incidence/RESEARCH_POLICY.md → NOTES.md → papers/poisson_subvariety_incidence_ja.tex。GitHub main唯一正本、TeX・コード・入力・証明書追跡、生成PDF非追跡、GitHub Actions禁止、Libraryは実push失敗時の緊急退避のみ。

**上位目標**：Poisson不変超曲面・完全交叉・葉層・退化五次因子・余法線Lie構造を通じた幾何学的構造定理によって、既存147被約成分の由来・境界・交叉を説明する。分類の完全性を定理の仮定に使わず、既知文献と本研究の貢献を分離する。

**証明済み**：smooth cubic三次Fanoの任意Poisson構造はrank2に一意延長。smooth degree≥4 hypersurfaceは非零ambient Poisson不変にならず、rank4退化五次因子の滑らかな既約成分次数は≤2。四次del Pezzo (2,2)の一般(S,σ_h)はrank2の被約一点延長、20変数Jacobi表示と28次元Rat122 incidence。smooth quadric一般Rat/Log/Aff延長(一点/A1/二点)と特殊Log境界。

**新しい構造定理 v18**：C=Z(σ_h)は反標準零曲線。rank4延長のPfaffian DはI_S I_Cに属し、法線一次ジェットはσ_h·κ（κは余法線Lie括弧、O_S(2)^2切断）。可換κ=0 ⇔ D∈I_S²（必要五次線形空間15次元）、非可換κ≠0ならDはSの一般点で一次消滅（従来必要空間41次元）。非可換κの消滅集合はCに含まれ、Cが滑らかな楕円四次ならC全体または長さ16の零次元集合。存在・十分性は未証明。

**具体的な一手**：既存141階数4源の零スキーム・PfaffianのS沿い重複度と余法線Lie括弧κを使い、滑らかな(2,2)曲面保存の候補を絞る。最初に(1)固定C楕円四次のA10 lift、(2)Pfaffian I_S²かI_S I_C＼I_S²か、(3)κの非可換化と16点条件、(4)20係数Jacobiの全解、を検査。一般点だけで特殊点を排除しない。次点は特殊smooth quadricの全scheme延長と既存147源の対応。

**重要な未解決**：特別なsmooth del Pezzoのrank4延長の存在／普遍非存在、特殊quadricの全分類、幾何構造から147を再構成する定理、一般次元Fano。副成果の完全等方Log点の非被約接錐・局所次元28・交叉23/24等は papers/isotropic_log_tangent_cone_ja.tex と旧現行計算証明書に保持し、削除しない。元のambient局所環の被約性は未証明。

手動再生成：bash scripts/verify_github_source.sh incidence（13 producer＋保存JSON＋LuaLaTeX）。組版と全計算は実行した分のみ報告する。最新の数学証明はpapers/del_pezzo_pfaffian_incidence.tex。研究一式に旧ZIP・作業ログ・PDF中間物は含めない。
