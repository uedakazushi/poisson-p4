# 固定T38フレームでの正weight6 edge ansatzの有界検査

2026-11-03。先に凍結した `Q_T38_QA8OVER9_COMPACT_FRAME` の四次体と実固有フレームを使う。
今回の結論は一つの明示的なedge valuation ansatzが成功しないという限定的なものである。
T38の全閉包から正軸やP2を排除したとは主張しない。

T38の四つのlinear smoothing tensorは、それぞれ

    x2*x4 d0 wedge d1,  x1*x3 d0 wedge d2,
    x4² d1 wedge d3,    x3² d2 wedge d4

で、物理weight `W=i+j-p-q` は `w=(-5,-2,-4,0)` である。
これらを固定したS0でpushforwardして、実canonical150×4係数行列Eを構成した。
各tensorのdivergenceは零のままであり、四本とも最高weightは7である。

S_ij(t)=S0_ij*t^(j-i)というフレーム則から、target rowのweightをWとすると、
pushforward係数は正確に `E(row,j)*t^(w_j-W)` となる。
正weight6項はW=-1なので、最初に試すaligned ansatzは

    tau_j=c_j*t^(-1-w_j),  c_j in K.

このlinear edge部分はrowごとに `sum_j E(row,j)c_j * t^(-1-W)` となる。
W>=0の全rowを取った行列のK上rankは4である。非零4次minorと逆元、QQへの非零normを保存した。
従ってregularな極限のために必要なW>=0のpole消去条件はc=0を強制する。

この必要条件はfull T38 sourceのcross項を加えても変わらない。
source periodはp0=mu*gamma12*t^5であり、各cross tensorのweightは対応する二つのwの和である。
従って同じ固定target rowで、toric部分、linear edge部分、cross部分のt冪はそれぞれ

    5-W,  -1-W,  -7-W

となり、互いに異なる。同じ係数の異なるt冪はpoleを相殺できない。
従ってこの固定frame・定数amplitudeのaligned ansatzでは非自明なregular曲線は得られない。

別のsource seed、固有frameのt依存摂動、非aligned valuation、複数段階のedge調整は今回扱わない。
特にP2の実枝、正軸のT38所属、他の24次元MPS成分の所属は依然未解決である。
新たな八変数高次jetや大規模GBは実行していない。

code: `code/q_T38_Qa8over9_edge_weight_test_20261103_v63.py`。
report/certificate: `reports/Q_T38_QA8OVER9_EDGE_WEIGHT_TEST_20261103_v63/`。
