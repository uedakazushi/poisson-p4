# 元の product41 の不良集合は有限

2026-11-03。標数0の開

`FULLU={a*P*F*B0*C*J!=0}`,

ただし `P=ad`、`F=500a+90e²+3e−24`、
`B0=(50a+18e²−9e−2)P+(60e−40)a`、
`C=100a+18e²−9e−2`、`J=5P−6e−4` とする。

**結論。** 元の六二次式に元の七 light 変数を掛けた42本の実際の多項式の階数は、
FULLU の有限集合を除く全幾何点で正確に41である。階数41未満の閉集合は零次元以下である。

証明は二つの独立受理済み有限入力の合成である。

1. `S61_NONSQUARE_PAIR_FINITE_BAD_LOCUS_20261103_v63_REVIEW` は actual804係数、
   pair chart、QQ Bezout 恒等式と先行 gcd/binary reduction の依存により、
   actual RAW24×56 の階数24未満の集合が FULLU 上で零次元以下であると証明する。
2. `S61_NONSQUARE_RAW24_PRODUCT41_BRIDGE_20261103_v63_REVIEW` は、
   元の物理六式の unit 行変換、実際の ell 座標と (Z,W) 商での四式の独立性、
   17 kernel 方向と非零 syzygy を独立に検証し、各幾何点で
   `RAW24 rank24 => original42 rank41` と証明する。

従って元の product rank<41 集合は RAW24 rank<24 集合に含まれる。
全集合の scheme 構造や長さを決定する議論ではない。有限集合が空、全 FULLU で rank41、
J=0、full cubic rank57、quartic 被覆、ambient 除去は主張しない。
新しい全開 pair-empty 証明が進行中であっても、この凍結証明書へ混入しない。
この corollary は構成的な pivot を新たに実行したという意味でもない。

検査入口 `code/s61_nonsquare_product41_finite_bad_locus_20261103_v63.py` は、両受理済み
report の正確な scope と入力・証拠 hash を照合する。数学的合成のみであり、
GB・新しい有限体・大行列計算は不要である。
