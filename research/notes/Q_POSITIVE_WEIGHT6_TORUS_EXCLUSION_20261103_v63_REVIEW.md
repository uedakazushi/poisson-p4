# 全 period・edge source の torus と閉 action-rank による排除の独立監査

2026-11-03。`Q_POSITIVE_WEIGHT6_TORUS_EXCLUSION_20261103_v63` を承認する。
正 weight6 軸の全非零点は actual T01、T02、T03 の全閉包の外にある。
T01 は別の全274式による排除とも一致し、T02/T03 が今回の追加排除である。

diagram の文字の数だけを数えず、一般 period を全て残した toric tensor を
線形 Pfaffian dual から独立に組み立てた。T02/T03 では唯一の edge の実二次 monomial を加え、
保存された全 symbolic raw source と全150係数で一致させた。どちらも Jacobi を満たす direct source であり、
未考慮の高次 smoothing correction を落としていない。

toric 項 xi*xj*partial_i wedge partial_j は全対角 torus で固定される。
二つの実 edge の character はそれぞれ

    (1,1,0,0,−2),    (1,1,0,−1,−1)

である。各 character の和は零なので PGL の四次元対角 torus へ降り、
kernel の連結 torus は三次元である。T01 では四次元全体が残る。
weight の代表として w0=0 を取り、その kernel の独立な全三方向、または T01 の四方向を作り、
全 symbolic source に対する有限 torus character の消失と同値な diagonal Lie action を
多項式微分で直接検査した。period や edge の一標本の stabilizer を一般化していない。

したがって source の全 admissible 点の projective orbit 次元は T01 で20以下、
T02/T03 で21以下である。PGL 変換で stabilizer の次元は変わらず、
既承認の full source の稠密性により、同じ上界は各成分の稠密部分で成り立つ。
canonical な線形 tensor 表現では、GL infinitesimal action と scalar 列を並べた行列の
rank 上界は小行列式の消失という閉条件である。その rank から一本の scalar 方向を除けば
projective orbit 次元になるため、この上界は成分の全境界へ延びる。
raw source が Euler wedge を含む場合でも、実 projective torus stabilizer の主張を使用し、
raw 150 行列の余分な垂直方向を projective rank と数えていない。

独立に承認済みの正軸の projective orbit は22次元であり、GL+scalar の23次 minor が非零である。
従って上の全閉成分の action-rank 条件と矛盾し、全非零軸は三成分の外にある。
T04以降や他の既知成分に同じ上界があるとは主張しない。

fresh producer と証明書も一致した。独立コードは
`code/q_positive_weight6_torus_review_20261103_v63.py`、receipt は
`reports/Q_POSITIVE_WEIGHT6_TORUS_EXCLUSION_20261103_v63_REVIEW.json`。
