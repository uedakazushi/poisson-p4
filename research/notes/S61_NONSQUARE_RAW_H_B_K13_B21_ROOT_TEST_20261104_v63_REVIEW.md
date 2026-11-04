# B21(d=s) と固定全35小行列式の独立監査

凍結済み B21_ROOT_TEST を、直前の独立 d-Euclid reviewer が持つ別実装の Ka13 固定
monomial-table 演算で再生する。作者の Ring、係数逆元探索、追加minor探索は使わない。

全84係数は独立受理済み source bridge の元の B21 係数へ結び戻す。全21 entry の
三つの Horner stages と、B_entry(d)=(d-s)Q+B_entry(s) の四係数恒等式を確認する。
各21値は13項、計273 QQ 項である。その7×3行列について、保存された35 row triples
すべての Laplace2×2 cofactors と determinant を別演算で検証し、すべて厳密に零となる。

体性・既約性・reducedness を仮定せず、Ka13全体での statement である。ただしこの
監査は prepared B の必要 candidate を保持する段階であり、original minor ideal と
FULLU open の明示的な scheme map を別に受理するまでは actual Bbad の存在・空性を
結論しない。full cubic57・quartic239にも移さない。結果と入力 hash は report に保存する。
