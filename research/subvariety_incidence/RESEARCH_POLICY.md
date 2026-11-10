# Poisson不変部分多様体と147成分分類の幾何

2026-11-10 / v6。再開順序：ルートRESEARCH_POLICY.md → REPRODUCE_GITHUB.md → 本方針 → NOTES.md → papers/poisson_subvariety_incidence_ja.tex。

## 価値最大化の原則

研究方針の上位原則は数学的重要性・新規性・概念的説明力・他成分への波及・正確性・実現可能性の総合的な最大化。新しい証拠や文献に応じて研究目標を柔軟に変更する。元目標の解決点・残る課題と変更理由を現行文書へ整理するが、履歴は蓄積しない。

## 主要結果と証明状態

- smooth Q3の一般Rat/Log/Aff intrinsic三型のPoisson延長ファイバーは、それぞれ被約一点（rank2）、被約A1（全rank2）、被約二点（rank2とrank4）。一般元の結論であり、特殊元の全分類ではない。
- 一般rank4のsmooth Q因子の重複度は、strict q-Casimir liftのmodular場がintrinsic葉に接するかどうかで判定できる。q²|D ⇔ v|Q∧σ=0。Q+C3と2Q+Hの幾何学的な違いを説明する。
- **新しい統一定理**：共鳴Log型 `P0=2ξ∧R+ξ(qU)J`（ξの重みは1,1,-2、q=qU+w0²+w1²）を考え、qUを非退化とする。重み1部分Hの階数は2か1に限られる。階数2なら、交差一次式b≠0の全延長は被約A1（rank2のみ）、b=0ならSpec C[t,v]/(v(t+2))（rank4枝）。階数1なら標準形qU=u0²+2u1z+cz²で、c≠0なら被約A1（rank2）、c=0ならSpec C[t,v]/(v(t-1))（rank4枝）。rank4の実退化因子は階数2型がQ+V(z)+V(cz²+w0²+w1²)、階数1型がQ+V(u0)+V(u0²+w0²+w1²)。この二種類の余次元二特殊層が非退化qUの全場合を尽くす。数学的証明は `papers/quadric_log_orthogonality.tex`（現行主稿からinclude）に保存し、零スキームの幾何とstrict stabilizer、全Jacobi理想を用いる。
- 一般4次del Pezzoの延長剛性、滑らかな三次の一意延長、rational成分の次元36/28、単純超平面の横断対称性の必要十分条件と既存147成分分類も保持。
- 計算は既存五実例とrank1の二正規形のJacobi/Groebner/Pfaffianを独立SymPyコードで厳密検査。標本計算を一般証明と同一視しない。新規性は一次資料による網羅調査が未完了であり、既知結果の再取得と本研究の貢献を区別する。

## 未解決の核心と優先順位

今回、**共鳴Log源のqU非退化域**をrank1/2 Hを含め全分類した。これは特殊Log全体、Aff特殊ファイバー、147成分との一対一対応の完全分類ではない。次はqU退化の境界やspecial Affとの交叉を研究し、rank4枝の所属と退化限界の機構を理解する。個別の複雑な場合分けが価値を下げると判断したら、4次del Pezzoの真の例外延長、P3横断Poisson対称性、147成分の概念的統一を優先する。

## 再現・保存

GitHub mainを唯一の正本とし、TeX原稿・共有TeX・7本の正確なSymPy producerと7つのJSON証明書・再現手順を保持。手動 `bash scripts/verify_github_source.sh incidence`。PDFはGit追跡しない、LuaLaTeX/latexmkで生成。GitHub Actions禁止。Libraryは実際にGitHub pushが失敗した場合の緊急退避に限り、通常保存しない。
