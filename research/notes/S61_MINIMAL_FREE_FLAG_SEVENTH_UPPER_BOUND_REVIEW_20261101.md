# free-M0 七次上界の独立証明・source監査

2026年11月1日。[主証明](S61_MINIMAL_FREE_FLAG_SEVENTH_UPPER_BOUND_20261101.md)
とproducerを読み、有限切断行列の非零minorから関数体の上界を得る論法を受理した。
独立のraw7再実装、二つのclosureの再実行、全Macaulay minorの保存・再検算は行っていない。
入力hashと次数・階数の算術を別checkerで確認した。

## 確認した範囲

元135式から定数項とJacobianを再構成し、保存したJとのQQ一致を検査している。
J T=0と選択J U=Iを確認してから、同じ元二次項から陰関数展開を七次まで生成する。
最高次も含め途中の初期idealへの縮約をしていない。
各有理係数のp=32003での分母可逆性を確認するので、有限個の分母を可逆化した
Z_(p)[a0,…,a4]上で同じ有限ジェット行列が定義される。
numpyの内積・bincountに対する整数上限も、int64および厳密な倍精度整数範囲内である。

二つの有限体計算は同じactual jet入力の切断idealを対象とする。
filtered closureは行演算と変数倍を閉じ、もう一方は局所standard basisの切断計算である。
主証明書の一致した階数27489は31824列の有限行列に非零minorを与える。
従ってQQ(a0,…,a4)上の同じ切断idealの階数も27489以上であり、
全商長は4335以下となる。有限体の一般fiber幾何を標数0へ移す議論は使わない。

既に証明された真六次までの一般Hilbert関数の総和は2568である。
したがってHF7≤4335−2568=1767となる。
既存J6の七次商1797は主producerの別QQ(a1) Gröbner計算に依存し、
これとの差から新しい七次方向が少なくとも30個存在する。

一般HF7が1767と等しいとは結論しない。
全初期ideal、実完成環、t冪零性、飽和除去へは進まない。

## 再現

```
.venv/bin/python code/s61_minimal_free_flag_seventh_upper_bound_review_20261101.py --output /tmp/free-flag-seventh-upper-review.json
```
