# 同じT38フレームでQ_(8/9)を保つedge seedは零だけ

2026-11-03。固定した四次体フレームで、別の非零constant seedを先に選ぶ方法を有界に検査した。
source weight w=(-5,-2,-4,0)に対し、tau_j=c_j*t^(5-w_j)と置く。
sourceのperiodはp0=mu*gamma12*t^5であり、全三つのcross correctionも同じlevel5となる。

pushforward後の物理係数は、各rowのweightをWとして、Qaに到る既存toric項に
`t^(5-W)` 倍の次の多項式を加えた形である。

    linear_edge_row(c) +
      [-2*cross03_row*c0*c3 + 2*cross12_row*c1*c2 +
        4*cross23_row*c2*c3]/(mu*gamma12).

regularかつt0で**同じ物理的Qa_(8/9)**となる必要十分条件は、W>=5の全8係数が零であること。
raw rowは(90,105,120,121,135,136,137,140)である。
八式は7 monomial `(c0,c1,c2,c3,c0c3,c1c2,c2c3)` の定数K線形結合である。
最初の7式の7×7係数行列はKの単元行列式を持つ。全逆行列と行列式・逆元・QQへの非零normを保存した。
従って各c_jは元の7式の**定数係数線形結合**であり、逆に全8式は(c0,c1,c2,c3)に入る。
これにより全8式のidealは正確に(c0,c1,c2,c3)である。凍結証明にGBは必要ない。

同じ固定frameで非零constant seedを先に選び、その周りのlevel−1 edgeへ進む経路は得られなかった。
正規化を変えるt依存群作用、別のperiod/eigenvalue frame、非aligned valuation、複数段階のsourceは扱わない。
T38の全閉包、正軸の所属、P2の実枝を排除したものではない。

code: `code/q_T38_Qa8over9_seed_weight_test_20261103_v63.py`。
report/certificate: `reports/Q_T38_QA8OVER9_SEED_WEIGHT_TEST_20261103_v63/`。
