# 補助五二次式の第4の安全な unit pivot（2026-11-03）

これは `Q0,Q1,Q2prime,Q3,D` の係数行列の rank 5 を目指す補助計算であり、主三次行列の constructive pivot 34 とは別である。

既存の安全な3 pivot 後の実際の QQ 係数行列は 2×18 である。第0行に対し、各係数から既知 unit `a,P,J=5P-6e-4` の因子だけを除いた後、係数次数2までの有限線形問題を解いた。modular rank は候補選択にのみ使い、受理は元の18係数と保存 cofactor の **全 QQ 多項式恒等式** `sum(A[0,j]*c[j])=a^3*P^3*J^2` の検証による。保存 cofactor は18個すべて非零、合計175項、元の行へ戻した後の最大次数4である。

この組合せを新しい第18列として append し、上の既知 unit を pivot にとった。新列は既存列の組合せなので image は変わらず、unit pivot の Schur 消去は localized cokernel を保つ。実際の native 出力は 1×18、8475項、次数19、最大分子/分母 bit 長510/528である。この第4 pivot は expanded 形式でも保存・再現される。未知因子は一切逆にしていない。

残った18係数はそれぞれ QQ 定数と既知 unit 因子だけで割り、すべての元係数との恒等式を検査した。primitive 行は7865項、次数18、最大541 bit の整数係数、1,291,970 bytes である。primitive cofactor 次数0,1,2、既知 unit monomial の因子数8以下という有限探索はすべて不成功だった。これはこの探索クラスに限る不成功であり、localized scalar ideal の非 unit 性も、全次数での不可能性も主張しない。3探索の合計は約2.1秒、重い solver は走っていない。

したがって証明済みなのは補助係数行列の rank が全 `a*P*J!=0` 上で少なくとも4であり、全 rank 5 は未解決である。pair gate も別に未解決である。raw 24×56 の全 QQ surjectivity、元の product rank 41、full cubic 57、ambient removal はこの packet からは従わない。

再現: report archive を repository root に展開し、既存 Cloud activation 後に `python code/s61_nonsquare_five_quadric_fourth_unit_20261103_v63.py`。有限探索を繰り返さず、元の175項恒等式、native Schur 出力の全 byte、全18個の unit normalization 恒等式を検証する。report に source/cofactor/helper と出力の SHA を固定した。
