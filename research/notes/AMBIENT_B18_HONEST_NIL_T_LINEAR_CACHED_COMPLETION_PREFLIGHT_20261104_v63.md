# 一次層の最終cache gate：未知proper/gauge二次値を要求しない

2026-11-04。全次数5×849作用素の最後の681 X列に対する、追加のsource感度検査である。
実二次値・一次関係のrankはまだ計算しない。

actual ResATのrow20は原点で全48成分が零だった。従ってResAT0*Bold2は最初から零。
ResAT2*Bold0も既に最終681 constant sourceで零と証明した。
最初の640列から必要なのは既存ResAT1*Bold1だけで、旧正則商の二次値は一つも要らない。

残るnative rowhにはu5とcore cross quadricsが出る。未知proper二次値に由来する
core quadratic積のscalar出力は、pure tau originとnu(0)=0により零である。
u5の未知even二次値五つについても、全681 final constant sourceで感度が零だった。
これにより、未知gauge2を解いたりu5_even2を新たに回収したりする必要もなくなった。
未知u5の値を零と主張するのではなく、正確な全source functional上の相殺である。

保持すべきmonic二次値だけは、元DAGの11節点と四つの保存済み57-source parameter
sA,tA,sB,tBからnu=t+s二乗/4として得る。CB own-F9からCA-F9へのfull quadratic
adapterを適用し、共有physical9二次polynomialがCA側と一致することをworkerで確認する。
新しいsource/gauge/IFT方程式は解かない。

whole preflightの111,057 QQ saved+live上限をそのまま使える。さらに19,912節点の
physical AST二次評価とu5回収を省けるため、この追加source gateは予算を増やさない。
条件付き承認された一回のcache値workerで全849一次列を評価できる。
ここでは全一次層、総長、nil生成元、三乗零性を主張しない。
