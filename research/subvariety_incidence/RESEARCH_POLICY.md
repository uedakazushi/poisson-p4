# Poisson不変部分多様体と147成分分類の幾何

2026-11-10 / v5。再開順序：ルートRESEARCH_POLICY.md → REPRODUCE_GITHUB.md → 本方針 → NOTES.md → papers/poisson_subvariety_incidence_ja.tex。

## 上位原則：数学研究の価値最大化

数学的重要性、新規性、概念的説明力、他成分への波及、正確性、実現可能性に応じて優先順位を変更し得る。進展で方針変更する際は、元の問題で証明したこと・未解決なことと変更理由を現行文書に整理する。歴史ログを蓄積しない。

## 現在の結果（論文に証明を保持）

1. Smooth Q3のintrinsic一般三型のambient延長：Rat(1,2)は被約一点（rank2）、Log(1,1,1)は被約A1（全てrank2）、Affは被約二点（rank2/rank4）。Ratはstrict stabilizer rank10、Logはrank9、Affは実モデルと既知一般軌道を使用。特殊点への一般化は未完。
2. Smooth Q3のrank4延長のmodular重複度定理：β=i_(P∧P/2)Ω=f dq（deg f=3）、v=div P∈so(q)、D=2qf、Q²|D と v|Q∧σ=0 は同値。従ってQ+C3と2Q+Hを幾何的に区別する。
3. **今回の一般化：共鳴Log境界の直交性判定**。V=U1⊕L⊕W（次元2,1,2）、q=H(u)+2z b(u)+cz²+w0²+w1²、HおよびqU非退化、bは任意の一次式、cは0も可。ξ=(u0,u1,-2z)、P0=2ξ∧R+ξ(qU)Jとすると、全Poisson延長は b≠0 で被約A1（rank2のみ）、b=0 で Spec C[t,v]/(v(t+2))（rank4の第二枝）となる。実rank4因子は Q+V(z)+V(cz²+w0²+w1²)。零スキームの直線+平面二次曲線からstrict so(q) stabilizerを直接求める通常の証明を papers/quadric_log_orthogonality.tex に収録。条件b=0は二つの一次条件で余次元2、以前の一変数wallはその断面。c=0もqUが非退化なら定理に含まれる。
4. 滑らかな三次超曲面のrank2一意延長、一般4次del Pezzoの被約一点剛性、Rat113/122の忘却次元36/28、単純超平面の横断対称性判定、Q上の明示的rank2/4例も保存する。
5. 有限実例の厳密検算はcode/subvariety_incidence/check_log_orthogonality.pyと対応JSONで五つの非退化二次形式（c=0を含む）のstrict stabilizer、Jacobi、q-Casimir、Pfaffianを確認した。一般定理は標本計算からの推測ではなく零スキームを使う通常の証明。数学的内容の文献上の新規性・優先権は未確認。

## 現時点の未解決の核心と次の作業

今回、従来の一次元wallを余次元二の不変な条件へ拡張して解決した。一方、H退化／qU退化の境界では零スキームによる同じ標準形が失われ得る。次はこの境界を系統的に層別し、rank4枝の有無を実Jacobiで判定する。また、特殊Aff型との交叉・147成分への正確な帰属を調べる。
ただし特別な重みの細分化が数学的価値を下げる場合は、4次del Pezzoの真の例外延長、P3の横断的対称性、147成分の統一定理へ優先順位を変更する。目標は新しい有限表の数ではなく、分類全体の幾何学的理解である。

## 正本と再現

GitHub mainが唯一の正本。原稿・共有TeX・6本のSymPy producer・6 JSON証明書・再現手順を全てmainに保持し、手動で bash scripts/verify_github_source.sh incidence を使う。生成PDFはGit追跡しない。LuaLaTeX/latexmkを使用、GitHub Actions禁止、Library保存は実push失敗時の緊急退避のみ。証明、局所組版、コード検査、main保存を区別して報告する。
