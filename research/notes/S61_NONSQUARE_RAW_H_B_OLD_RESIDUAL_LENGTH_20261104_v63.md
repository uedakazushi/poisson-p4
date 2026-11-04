# 旧有限残余の正確なlength117とモニックH-cover length468

R=QQ[a]、A=R[e]/(h1)とする。h1はeモニック16次なのでAはR-free16。
保存したMはA上のk1乗算写像であり、そのcokernelはR[e]/(h1,k1)そのもの。
整数行列N=DMのDはQQ単元で、detNはa次数148、a-adic valuation31と厳密に
確認されている。

PID RでのSmith標準形を取ると、cokernelのQQ次元は各非零invariant factorの
次数和=deg(det)=148。aを逆にするとa-primary部分だけが消え、そのlengthは
v_a(det)=31。従ってQQ[a,a⁻¹,e]/(h1,k1)は正確にQQ-length117。
未知因子のfactorizationやreduced性は不要である。

この有限代数上でも既存Hはdモニック4次なので、old containing cover
QQ[a,a⁻¹,e,d]/(h1,k1,H)のQQ-lengthは正確に468。
これは実際のB-bad schemeそのもののlengthではないが、旧残余の幾何学的点数を
高々468に抑える。LとQの全線形成分は別のunit証明で除去済み。残るQ因子の次数和32
による高々256点と合わせ、全B-badの幾何学的点数は高々724。

724はscheme lengthの主張ではない。他のcover議論はsquarefree化と零点集合の包含を
使うためである。117/468のexact lengthは上記明示した有限商代数だけに適用する。
新たな点・素数・GB・solverは使わない。
