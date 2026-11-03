# 正 weight6 軸の projective 軌道と追加排除の独立監査

2026-11-03。`Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63` を承認する。
全 eta≠0 点は一つの22次元 PGL5 軌道であり、T39、T40 と、既知の22次元二成分
H+2Q_smooth、3H+Q_rank4 のいずれにも属さない。他の成分や全局所枝の尽くしは結論しない。

実 tensor P_eta の成分を独立に組み立て、eta=1で各 E_ij に対する微分を

    E_ij*P+P*E_ji−x_j*partial_i(P)

から直接多項式微分した。全25列・150行と scalar 列が保存行列に完全に一致した。
全方向の divergence も零であり、canonical tensor の線形部分空間にある。
Euler wedge の quotient による見掛けの方向を rank に含めていない。
23次非零 minor と、全行列の三つの独立な kernel ベクトルを厳密に検証した。
後者の scalar 成分はすべて零なので GL25列の rank は22以下である。
scalar を加えた23 minorからGL列の rank は22以上となり、正確に22である。
scalarを加えた rank23から一本の affine scaling を除くと、projective 軌道次元は22になる。
GL rankだけを projective rank とみなす推測ではない。

有限変換も全150係数で直接検査した。
g(t)=diag(1,t,t²,t³,t⁴) の pushforward を t5 で割ると P_(eta/t6) に厳密に一致する。
複素数上で任意の二つの非零 eta の比は六乗であるため、非零軸は一軌道である。
実際の多項式曲線 eta→P_eta が eta0へ延びるので、この軌道の閉包に Q_(8/9) が入る。
既承認の全域 graph 判定では Q_(8/9) は T39/T40 の外にある。
どちらも PGL 不変閉成分なので、軌道の一点が含まれればその閉包の極限も含まれて矛盾する。

二つの quadric 成分の実成分性・次元22は、凍結済み v03 と v05 の証明を依存として再利用した。
もし一方に22次元軌道が含まれれば、その閉既約22次元軌道閉包は成分全体になる。
一方、五次退化式を与える polynomial map の下で、この軌道は全点が純五乗に写る。
純五乗は Veronese の閉錐であり零も含むため、その逆像は閉で、軌道閉包全体を含む。
既知二成分にはそれぞれ H+2Q_smooth、3H+Q_rank4 の非純五乗を持つ実一般元があるので矛盾する。
退化式が零になる境界で projective な退化式写像が未定義になる問題は、閉錐の逆像を使うことで避けている。

fresh producer と全証明書も一致した。独立コードは
`code/q_positive_weight6_orbit_review_20261103_v63.py`、receipt は
`reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW.json`。
今回の結論と、別途承認済みの T01/CA/CB/N 排除は区別して保存する。
