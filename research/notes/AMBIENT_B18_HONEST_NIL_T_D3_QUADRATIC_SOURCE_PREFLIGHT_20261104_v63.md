# d3だけの153二次source：十二adapter callbackへの感度

2026-11-04。全一次rank35から得た既知二次rank174の六商へ、d3だけのsourceを射影する。
実d3一次行列のrankは15なので、全constant kernel153本を固定minorで構成・照合する。
各sourceは真に一次零であり、その二次値は初期関係となる。値自体はまだ求めない。

六商dualはwedge三targetをr二乗でしか読まない。従って
−C1 Wunit0^-1 d3_pivot,1のwedge部分に必要なのは、81pivotのC1のr係数だけとなる。
これはrow20 X2のz0r,z4r,z5rで定まり、保存済みB11に属する方向である。
top二targetのC1は既存X1で全方向が定まる。

ResAT0=0なので旧Bold2は不要。残り得る未知は四zero graph64,68,73,74のadapter二次値を
三方向zi*rで読んだ十二callbackだけである。その係数を独立変数として保持し、
全153source・六dualの感度を固定QQ行列から計算する。零なら別のcache値workerへ進めるが、
非零ならliteral B11 proper sourceからその四係数を評価するsource gateが必要となる。
この検査で未知callback自体を零に設定しない。

X3、m2、zi*r二乗の新しいsource/gauge jetは、このd3限定gateには不要である。
元の814 source全体や真のHF2を決めるという主張はしない。

固定QQ検査は0.24秒で完了し、全153 source・六dualで十二callbackへの感度は全て零だった。
従って追加のsource/gauge値を求めず、既存B11のnative値とA1/B1・X1・m1だけによる
別の小さい値workerが可能である。ここではその値や追加rankをまだ計算していない。
