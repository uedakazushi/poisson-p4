# K92 の非零 determinant を元 binary source と別の係数計算で確認

2026-11-05。主計算の整数多項式 Bareiss を再利用せず、凍結済み K92SP001 の binary 仕様を独立に読み取った。完成済み source の literal 多項式だけを読み、保存されていた未完了 determinant DP は実行しない。producer の parser や計算関数も import しない。

元11対11行列について全 source の e,d 指数が非負であることを確認し、e=d=0 の係数面を取り出した。行 index 1 に index 2 の二倍を加える determinant 一の操作も独立に施した。保存された terminal row/column potential と各行の整数分母から正規化し、全121成分・579整数係数が軸 source と一致した。

    Caxis = diag(L_i*a^r_i) Mprime(a,0,0) diag(a^c_j),
    sum r_i + sum c_j = −5,
    det M(a,0,0) = a^5 det Caxis / product L_i.

符号付き計算の前に、絶対係数の subset DP で次数7までの support・高さを固定した。全 owned 上界37296、係数548 bit、積101320対、入力4.23MBと新出力上界2MBである。実計算では47592対、470 bitで、

    [a^12] det M(a,0,0)
      = 488042709959022324503443014699328840826617856
        / 17514302730560302734375 ≠ 0

を得た。11行の permutation subset determinant であり、主計算の多項式除算や Bareiss patch は使わない。主計算の全24係数一致を検証したという主張ではない。

これは元の generic 多項式 P_K92=det M(a,e,d) の非零を証明する。d=0 は元 FULLU の局所化の外にあるため、その点で Poisson 行列が full rank という主張ではない。多項式の係数面による非零証明であり、新しい数値 beta 点を選んだものでもない。

旧 Omega と D(P_K92) の積は非零なので、この新しい generic chart は非空である。PIVOT14 の十分条件付き source 操作は残余1対8へ進められるが、全 FULLU57・quartic239・補集合の尽くし・実 bad Fitting ideal の同定は別課題である。

再現は [独立コード](../../code/s61_nonsquare_K92_axis_a12_peer_20261105_v63.py) を `--mode source`、次に `--mode value` とし、それぞれ新しい空の出力先を指定する。
