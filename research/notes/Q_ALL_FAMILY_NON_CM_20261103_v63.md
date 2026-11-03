# 全Q_aでの三つの既知成分と全Poisson局所環の非CM性

2026年11月3日。全有限複素aについてQ_aはT01、CA、CBに入る。
この三成分の射影次元は25、28、28である。従って全Poissonスキームとその被約空間の
Q_aでの完成局所環は、全aでCohen–Macaulayではない。

T01への帰属は別の明示的な可換行列退化
[Q_T01_MEMBERSHIP](Q_T01_MEMBERSHIP_20261103_v63.md)と独立監査に依存する。
CA/CBへの帰属を、rank4開からa=0へ不明瞭に延長する必要はない。
全proper模型に直接次の値を入れる。

源の座標を(h,x,y,z,w)=(x0,x1,x4,x2,x3)とし、

    H=(4/3)h*x*w−(1/3)h*z²−x²*z

とする。CAにはf=a*h²、g=−2/9を入れる。CBにはF=h²、
associated P1の比[a:−2/9]を入れる。この比はa=0を含め常に定義される。
fのx微分は零なので、二つのraw tensorは厳密に同じであり、既にdivergence零である。
上の座標順で戻すと全150係数が物理的Q_aと一致する。

Poisson空間の既知の実際の既約成分T01とCAの次元は異なる。
閉点での各成分の局所次元はそれぞれ25、28であり、完成後もその次元を持つminimal成分が残る。
完全Noether局所Cohen–Macaulay環はequidimensionalなので、全完成局所環はCMではない。
同じ議論は全Poisson空間を被約化してから完成した環にも適用する。
これは非被約性を前提にする議論ではない。

全局所枝の数、全局所次元の上界、T01自体の正規性・枝数は未決のままである。
個々のCA/CBが全aで非CMであるとも、この証明からは結論しない。
従来の指定Q1での個々の成分の非CM定理と、ここでの全Poisson環の定理は区別する。

全係数確認と射程は `code/q_all_family_noncm_20261103_v63.py`、
`reports/Q_ALL_FAMILY_NON_CM_20261103_v63/report.json` に保存する。
