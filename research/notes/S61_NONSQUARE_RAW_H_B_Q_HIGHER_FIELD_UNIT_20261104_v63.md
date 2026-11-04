# 残るQ高次例外成分の逐次numberfield検査

factorindex6,7,8をそれぞれ別budgetで、K=QQ[e]/f（次数6,9,11）の8次元代数
K[a,d]/(Q,H)として調べる。元sourceは全22minorを保持する。
全q0/pivot逆元はQQ[e]のliteral Bezout=1で確認し、非単元なら停止する。

source/modelの全係数とserialized checkpointをsolve前に保存する。QQ全次元8rの
大きな行列は作らず、8×176 overKだけを使う。最後のcofactorの全係数identityを
独立quotient arithmeticで確認する。modelデータはcheckpointだけに保持し、native
mainに重複M/RREFを書かず、aggregate32MiBの実byte上限を維持する。

成功時のΣCiPi=U8により、その閉成分のFULLU部分は空。全共役e根・全a/d rootsを
含む証明で、point samplingではない。他成分への自動jobや、旧有限a残余の除去を
意味しない。点・新素数・GB・resultantは使わない。
