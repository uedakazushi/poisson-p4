# 全99関係の独立Macaulay2 elimination検証

2026-11-03。CloudのMacaulay2 **1.24.11** で、新しいformal normal modelの99関係が全kernelであることを、source変数tの全eliminationと比較した。
既存のセットアップは再実行していない。QQの`Algorithm=>LinearAlgebra`は使わず、default Groebner/eliminationを使用した。

F=nu−s(s−r)を新しいbase座標とし、odd idealの生成元を

    z1,...,z6,F,wr,ws,sr²

と取った。root側では二つのcore syzygyをs*Y7−r*Y8、r²*Y8−w*Y9で選び直した。
作者packetの第二syzygyは、これに第一syzygyのr倍を加えたものに等しい。他の43 syzygyと55積関係は、変数名とF置換後に全て一致する。

44 linear-in-Y syzygyと55積関係をKとし、source ideal

    (t²−(F+s(s−r)), Y_i−b_i*t)

からtを消去した全idealをEとした。M2で **K==E** が成立し、trim後の生成元も99個と出た。
これは全kernel比較であり、各関係のpullback零だけを確認した検査ではない。
全polynomial idealの等号なので、Cへの体拡大と原点でのcompletionにも移る。

実CA/CBをこのnormal modelへ移す形式幾何の証明は別の `AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63` に依存する。本検証だけでその適用を証明したとは扱わない。
M2 producerは60秒・約2.2GBの上限内で正常終了した。時間の実測値は保存していない。

再現は別の出力先へstdoutを保存する。

    M2 --script reports/AMBIENT_ORIGIN_FORMAL99_M2_20261103_v63/verify.m2

成功出力は `FULL99_EQUALS_ENTIRE_SOURCE_ELIMINATION_PASS` と `99` である。
