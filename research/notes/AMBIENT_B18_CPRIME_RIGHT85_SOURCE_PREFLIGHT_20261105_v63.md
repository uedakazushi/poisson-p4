# Cprime の一次 syzygy 85 本からの既存二次値のみの三次検査

V0=0、V1 の階数は 24 である。固定 24 独立な元列について、全九変数の一次倍は 216 列を持ち、既知二次 initial 空間の階数は 131 である。従ってその定数係数 kernel は 85 次元である。本 source builder は新しい大きい逆行列を作らず、疎な pivot と元 216 列の係数表示、全 85 kernel の自由行 I85 と零 multiplyback を保存する。

85 個の homogeneous な一次係数ベクトル u1 は、実際の元 V の選んだ 24 列に直接掛かる。V1 u1=0 なので、実関係 V u1 の三次 initial は V2 u1 となる。定数 kernel 列の三次計算と違って、V3、二次 cancellation lift、変数依存の列変換の高階係数は一切現れない。

既知 rank 26 の残余三次汎関数は pure-plane r,ell のみを読む。そのため u1 の z 方向係数は、未知の V2 の値にかかわらず寄与しない。r,ell 方向の係数に対しても V2 の rr,rell,ellell 三成分だけを読めばよい。これらは完全二次値 packet の全 44 列について既に確定している。ここには fresh source3、gauge3、未保存 z² 成分の零埋めはない。

この packet は全 85 adjoint と source kernel のみを保存し、既存の二次値との pairing はまだ行わない。いずれかが非零なら、実際の三次 initial が既知 rank 26 を 27 に上げ、Cprime/ell³ の長さ上界 9 と既存下界が一致する。全零の場合はこの 85 本による追加が無いことだけを意味し、残る定数 kernel 19 列や実 exponent の結論は出さない。

Source と将来の値 job は 60 秒、2 GiB、60000 owned numeric slots、2048 bit、入力と出力の合計 32 MiB 以内とする。値 worker は別の frozen code/custody と独立 SOURCE peer、親の GO の後にのみ走らせる。
