# Poisson部分多様体incidence：現行研究ノート

2026-11-10 / v2。証明の全文は `papers/poisson_subvariety_incidence_ja.tex`。研究再開の入口は `RESEARCH_POLICY.md`。

## Quadric上のmodular場の意味

非退化二次qのstrict Casimir lift `P` に対し `β=i_(P∧P/2)Ω=f dq`、deg f=3、`D=2qf`。`v=div P` には `v∈so(q)` と `[v,P]=0` が成立する。`dβ=-i_(v∧P)Ω` をEuler収縮すると、Q上の `v_Q∧σ=0 ⇔ q|f` を得る。rank4ならQの重複度は1または2。v=0ならdf∧dq=0からf=0、rank≤2。一般論は論文 `thm:quadric-modular-multiplicity` を参照。

**一般Ratの新しい厳密な開条件**：q=Σxi²、L=x0、
`R=3x0²+2x1²+5x2²+7x3²+11x4²+13x1x2+17x2x3+19x3x4+23x1x4`。
J(q,L,R)を保存する `so(q)` の10変数線形行列はrank10、nonzero minorは27980515917168640。既約Rat(1,2)源上で開条件なので、一般元のstrict stabilizerは0。本文の核 `K_T` と `δK_T=-3/2 v_T` により延長はJacobian一点。
exact producer `code/subvariety_incidence/check_quadric_symmetry.py` はこの行列と、`Q+C3` の `f=-t u³/4`、`2Q+H` の `f=-5x0 q/4` を検算する。標本は普遍定理の証明ではない。

**今後の障害**：strict automorphismが非零でもJacobi条件を満たす延長が存在するとは限らない。Log、Affの一般点の `aut(Q,σ)` を求めた後、その部分空間に10変数Poisson Jacobiを制限し、非被約fiber・特殊parameter・閉包を検証する。既存Log一般A1、Aff一般二点の結果は元の証明の仮定を確認してから使用する。これが高価ならdel Pezzo例外集合へ転換する。

## Del Pezzo例外層

CI22では全延長は `J_h+q0 A+q1 B` の20係数。一般pairについて普遍二次rank60、90 syzygy、二nonzero minor `70875/88`、`4443337729552/7688671875` から被約一点。特殊条件におけるminor消失は、真の延長増加の十分条件ではない。反標準hとpencilの五つの退化memberの位置を用いて層別する。

## 全分類との接点と一次文献

Cubic-incidence 40→Rat113次元36（忘却fiber4）、CI22-incidence 30→Rat122次元28（忘却fiber2）の証明は現行論文にある。P3の単純H判定は横断的strict対称性が必要十分、重複Hは別。LPT published DOI 10.1002/mana.201100354とMayanskiy DOI 10.1016/j.jalgebra.2014.08.049を一次資料として使用し、今回の新主張の文献上の優先権は未確認。

GitHub mainのみを正本とする。GitHub Actionsは禁止、Libraryは実push失敗時の緊急退避のみ。
