# B18 nil 三乗の有限原点帰着：独立監査

2026-11-04。MS の同名 primary の全依存 hash と証明を独立に監査した。
対象は実 (D=CA\cap CB\) の一般 B18 点の完成環であり、全 Jacobi scheme ではない。
新しい係数、source jet、点、行列 rank を計算していない。

## 完成と homogeneous 成分

relative weighted model の零次数係数環は既に完備であり、正次数変数の個数は有限、
全 ideal は Noetherian 完成環で閉じている。
nilradical は continuous torus に安定である。
重み (0,\ldots,N-1\) の finite Vandermonde projection を torus translate の線形結合で作ると、
任意の nil 元の一つの homogeneous 成分を正次数 (N\) 以上の tail を除いて得る。
その tail は maximal ideal の (N\) 乗へ入るため、閉じた nilradical に極限成分も属する。
したがって primary の homogeneous nil 成分を使う議論は完成環でも正当である。

安定 tail の nilradical は全体として square-zero であり、被約 tail とは置き換えない。
正次数 (A_n\to M\) の injectivity は、全87 graph を持つ同じ actual model の受理済み結論である。
nil homogeneous (a_i,b_j\) の両像は nil なので、(i+j>0\) の積は tail で零、従って source でも零。
これは零次数と正次数の混合積も含む。有限係数和から無限級数へ渡すと

\[
 \operatorname{Nil}(A)^2=\operatorname{Nil}(A_0)^2,\qquad
 \operatorname{Nil}(A)^3=\operatorname{Nil}(A_0)^3.
\]

正次数の任意の元 (c\) に対して (b_0c\) は正次数 nil 元であるため、
上の二乗・三乗は正次数 ideal に annihilate される。
一般比較 kernel の正次数消滅を全域へ外挿する仮定は使わない。

Laurent tail の同定は有限 polynomial model を先に局所化する。
完成した (w=0\) 環へ (w=1\) を代入していないことを別 scope receipt で確認した。

## 比較 kernel と有限長

正規化の actual tensor 比較は同じ physical (R_0\) 上にある。
各正規化は integral・Spec-surjective であり、その性質は二つの base change 後も保たれる。
従って比較 kernel (K\) は nil ideal である。
whole normalized (E_4\) の nil 三乗零から、零次数 nil 三乗は (K\) に入る。

conductor 積は (K\) を殺し、四つの実 conductor curve は Noether parameter
\(\ell\) 上有限である。CA-critical と Pa は whole nil square-zero、
Qb は whole nil cube-zero、CB-critical は actual CA の外である。
ここで使うのは full local ring の式であり、非被約 conductor の被約点だけを調べた値ではない。
よって nil 三乗の (K\) 内の像は全 generic support で零であり、有限 DVR torsion になる。
これは相対原点に有限長で支持される。

## 正確な条件付き結論

最初に疑った rank の混同はない。新しい critical-contact packet は

\[
 \operatorname{rank}_D K=1+1+0+2=4
\]

を同じ actual 比較で別に証明し、source と代数の独立 gate を持つ。
これは別の \(\operatorname{rank}_D C=4\) から推測したものではない。
既受理の \(\mu_DK=\operatorname{length}T-4\) により、
\(\operatorname{length}T=8\) なら (K\) は DVR-free4 となり nil 三乗は零である。

actual Qb family の一般点では whole nil 二乗が非零で、family は (P_0\) へ特殊化する。
coherent nil-square ideal の支持は閉であり、faithful completion と群 product により失われない。
従って (P_0\) では nil 指数は少なくとも3、上の長さ8条件下では正確に3となる。
逆に三乗が非零なら深さ19だが、三乗零から深さ20を結論しない。

共通19群方向を戻すと、nil 三乗は有限長横断 module の
\(\mathbb Q[[G_{19}]]\) への completed tensor であり、非零なら有限自由。
その rank は未計算で、formal support の次元上界は cone19・projective18である。
新しい global associated closure、成分数、原点の深さ、三乗零そのものは主張しない。
