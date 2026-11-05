# 実cache値job前の係数高bound

既知入力の各source classについて共通分母Dと、全係数cに対する
|D c|<=Uを整数LCM/絶対値で求める。積は分母・分子上界を掛け、各出力の最大項数を
掛ける。和は分母LCMへ揃えて正の分子上界を足す。これらは値の符号相殺に依らない。

raw33積、nu=t+s_trace²/4、primitive V*Wv、bar-rho dot、全A23、実native ell first、
B30 inverseによるa1、全pivot補正とnormalized graphまで伝播した。
最大の分母・分子上界は538bitで、別jobの許可上限2048bit以内。
Phi/A23/不変量/graphの符号付き和の実値はこの証明で計算していない。
