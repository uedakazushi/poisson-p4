# 五つの既存 map origin からの actual quartic dual gate

既知 L35/Q3 の次数 3 quotient は三つの literal r³ classes（minimal head 0,3,4）で生成される。保存済み q3 は全 600 三次単項式をこの三個へ送る構成的 reduction で、差は既知の実 initial 関係に入る。

次数 4 の各単項式を「三次 divisor×最後の変数」に分け、q3 を適用する。次に x*(head_f r³)=r*(head_f r²x) として、内側の三次式へもう一度 q3 を適用する。二つの差はいずれも既知三次関係の一次倍である。従って全 actual gr4 は三つの pure r⁴ classes により生成され、この二段 reduction p4 は明示的な source witness を持つ。

これは p4 が既知 quotient 上で勝手に well-defined と仮定するものではない。actual gr4 への surjection を構成した後で、既証の actual dimension 1 を使う。

old Told は m⁴Told=0 と確定している。actual S-linear natural map を五 old minimal head の r⁴ 倍へ適用すると、degree 4 の initial は map origin M0 を r⁴ に掛けたものになる。高い map 係数は degree≥5 なので関与しない。p4 のもとで、この五つの actual 関係は保存済み二行列の 3×5 積として得られる。

もしこの小行列の階数が 2 なら、actual gr4 が一次元であることから、その kernel はこの五関係の span と一致する。三座標の左 annihilator と p4 の合成は真正な actual quartic dual である。階数が 1 以下なら別の source 関係が必要であり、階数 3 なら既知 actual HF4 と矛盾するため停止する。どちらも都合のよい階数に置かない。

本 packet は p4 と source witnesses を固定するだけで、3×5 の値・階数・actual dual・nil cube はまだ計算しない。新 source、gauge、jet は不要である。次の bounded rank job の code/custody と SOURCE peer、GO は別に固定する。
