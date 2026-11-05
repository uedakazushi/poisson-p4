# 十二角軸Bareissのbinary witness予算

2026-11-05。旧M12_AXIS_SOURCEの全JSON witness予測40.33MBは正しいFAILとして保持する。
新しい係数・axis determinantは計算せず、同じ固定計算のpublicationだけを変える。

入力degree9、minor degree108、minor1096bit、long-division等のtemporary2202bit、
全owned33275を保つ。全update係数20504 slotと十二pivotを、degree/sign/byte-length/
integer-magnitudeという有界binary recordで保存する。全Bareiss exactdivision identityと
row/column permutationを復元可能にし、quotient witnessを省略しない。
matrix/pending product/partial longdivisionとcursorのguard stateも同じ有界recordを使う。

既存十一角の固定pivot permutationを十二角へ拡張して用いる。
左上十一角は既存CaxisのQQ整数row単位倍なので、その固定pivotは同じsource上で非零。
最後のscalarが零ならaxis rank十一だけを言い、generic determinant零にはしない。

binary main7,899,784B、全guard state10,465,401B、入力3,874,058B、合計22,239,243B。
旧JSON FAILを上書きせず、32MiB内に全witnessと途中stateを保持できる予算である。
最終determinantだけなら15.7KB以下。整数productのsource上界は2,068,638回で、
時間完了を保証せず45秒soft/60秒hardを想定する。2202bitは旧2048bitより大きいため、
実workerの許可はこの明示boundに対して別に必要である。scalar workerは未実行。
