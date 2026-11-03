# 修正 shear と28本の多項式 kernel による短い rank98 証明の独立監査

2026-11-03。`Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63` を承認する。
generic な全 canonical Jacobi rank98を、有理函数の長い行消去を使わず証明できる。
一般二次 obstruction の rankや exact Artin長はこの packet から結論しない。

密性に必要な有限変換を最初に独立検査した。
元の canonical150全族を斉次 tensor に組み、

    h'=h, x'=r*x, y'=r²*y,
    z'=r³*z+(r³−epsilon)*h/(6*n0), w'=epsilon*r*w

の真の pushforward と scalar epsilon*r^(−2) を計算した。
epsilon=±1の両方について、全150係数がパラメータ変換
(n0,n1,n2)→(epsilon*n0,epsilon*r*n1,n2/r) の tensor と恒等的に一致する。
determinant は epsilon*r7、必要な開は n0*r≠0である。
したがって D(n0*n1) で epsilon=1,r=1/n1 を取る密な n1=1表示が正当である。
定数 yz bracket を変えてしまう純粋な diagonal 特殊化は使用していない。

この二変数全族の26本の実接方向を別に組み立てた。
全GL25方向と scalar、n0,n2 の実族微分を取り、指定された26列を選ぶ。
これらに二つの保存された一次 polynomial vector を加えた126×28行列が、保存 kernel と完全に一致する。
独立の物理350二次係数辞書から QQ[n0,n2] の全350×126 Jacobianを作り、
全350×28積が恒等的に零であることを別の多項式演算で検査した。
28次非零 minor と98次非零 Jacobian minorを正確な QQ で照合したため、
generic kernel は正確に28次元、Jacobian rankは正確に98である。
二つの追加方向は実26次元方向と独立で、generic 二法線 quotient の基底となる。

さらに元の350物理二次係数を二つの150次元方向に直接代入し、
u²,uv,v²の全350×3 Hessian列も独立に一致させた。
このデータの cokernel rank が1であることまでは、まだ証明していない。
二つの法線が存在するだけで一般 Artin長が5になるとの推論はしない。

短い producer の fresh 出力・証明書も凍結物と一致した。
独立コードは `code/q_new_family_polynomial_kernel_review_20261103_v63.py`、receipt は
`reports/Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63_REVIEW.json`。
先行する三変数・長い native rank証明とその独立 review は変更せず、別の短い exact 証明として保存する。
