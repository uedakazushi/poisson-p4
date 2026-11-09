# Poisson不変部分多様体による147成分の幾何学的理解

2026-11-10 / v4・現行方針。再開時はルート`RESEARCH_POLICY.md`、`REPRODUCE_GITHUB.md`、本書、`NOTES.md`、`papers/poisson_subvariety_incidence_ja.tex`を読む。

## 研究全体の価値最大化

数学的重要性、新規性、概念的な説明力、他の成分への波及、証明の正確性、実現可能性を総合して研究価値を最大化する。優先順位は固定ノルマではない。真の障害・文献・新定理に応じて現在の有効方針を書き直し、元の目標の達成点と残る課題だけ記し、作業履歴を蓄積しない。

## 主要な到達点と証明状態

1. **smooth Q3のmodular因子重複度定理**：strict q-Casimir lift Pのβ=i_(P∧P/2)Ω=f dq、deg f=3、v=div P∈so(q)。一般rank4ではD=2qfのq重複度は1か2で、q²|D ⇔ v|Q∧σ=0。本文に通常の一般証明を収録。
2. **smooth Q3の一般intrinsic三成分の全延長**：Rat(1,2)は被約一点（rank2）、Log(1,1,1)は被約A1（全点rank2）、Affは被約二点（rank2/4）。Ratのstrict stabilizer消滅はrank10 witnessと開条件、Logは普遍Jacobi/Pfaffianとrank9 witness、Affは既存の明示例＋LPT一般軌道を用いる。全特殊点へ拡張しない。主稿は`papers/poisson_subvariety_incidence_ja.tex`、Logの全証明は`papers/quadric_log_extension.tex`。
3. **特殊Log境界の交叉二直線**：一点の全延長スキームはSpec C[t,v]/(v(t+2))。rank4枝の退化因子はQ+V(x2)+V(x2²+x3²+x4²)。対数一形式は因子的零点を持つ。全Jacobi理想を計算しており、単なる接空間観察ではない。
4. **今回の新しい一変数定理**：q_c=Σxi²+2c x0x2 (c²≠1) と共鳴Log境界σ_cでは、c≠0の全延長は被約A1、rank2。c=0では上の交叉二直線とrank4枝。strict so(q_c)対称性行列の九次minorは13824c/(c²−1)^4であり、直交性wallでrankが9→8に跳ぶ。したがって**因子的零点や重みの共鳴だけではrank4延長を保証せず、quadricの直交性が本質的**。一般論は本文、正確な記号検算は`check_log_wall.py`。
5. **保持する先行成果**：滑らかなcubicの一意rank2延長、一般CI(2,2) del Pezzoの被約一点延長、rational成分の36/28次元、超平面の横断Poisson対称性。原稿・証明を削らない。
6. **文献**：Loray–Pereira–Touzet, Math. Nachr. 286 (2013), 921–940, DOI 10.1002/mana.201100354はQ3のintrinsic Rat/Log/Affを分類。Mayanskiy, J. Algebra 424 (2015), 21–45はsmooth cubic。今回のambient拡張に関する文献上の新規性・優先権は未確認。通常証明、計算による同一式検査、独立査読は区別する。

## 現在の最重要な未解決問題

**Q3の特殊intrinsic構造の全延長ファイバーの層別**。今回の一変数wallを一般のLog共鳴・zero-factor・quadricの直交配置へ拡張する。まずstrict `aut(Q,σ)` をLie代数の重み分解から計算し、その小さな空間上のJacobi理想を求め、ファイバーの既約性・非被約性・rank2/4・実Pfaffian因子・147成分への帰属を分ける。現状のwallの全族への必要十分性は**未証明**。

次はAffの特殊点との境界や、del Pezzo例外的延長・P3の横断対称性を比較する。147成分の概念的構造定理が価値の高い最終目標。Feigin–Odesskii型、特異cubic、一般次元も進展により順位を再評価する。

## 再現・保存

GitHub mainが唯一の正本。現行原稿・本方針・ノート・5つのexact SymPy producerと5つのJSON証明書・手動再生成手順をmainに収録する。再現は`bash scripts/verify_github_source.sh incidence`。GitHub Actionsは使わない。Libraryは実際のGitHub push失敗時の緊急退避のみ、通常保存しない。不要な旧原稿・作業ログを積み重ねず、重要な結果と必要証明を維持する。
