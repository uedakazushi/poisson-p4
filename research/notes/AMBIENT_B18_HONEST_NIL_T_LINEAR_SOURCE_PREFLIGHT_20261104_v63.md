# T_nil一次層のsource感度：一行・十五monomialへの縮約

2026-11-04。本packetは既存値に対するsource感度のpreflightで、実二次係数や全一次関係を
まだ計算しない。未知m1やadapter2を零に設定しない。

固定minimal frameの三wedge射影は、正確にsector0の(0,7),(4,7),(5,7)の座標選択である。
さらにm0の第7行は、元free44 source20における5/2だけが非零である。
従って三wedge一次値に必要なX2はrow20だけとなる。canonical mの第7行は
Vgamma_A(z=0,r)/rなので、m1はr方向だけを持ち、その値は保存済みC44のpure-r二乗係数。
第三sectorの未知pure-r三乗や全m1を計算する必要はない。

ordered Hadamardのprefixは、target iについてzi*zjのj<=iとzi*rだけを残す。
i=0,4,5の三targetでは合計十五monomialであり、そのうち三つzi*rは既存B11にある。
残る十二はz_i z_jで、現B11から任意に零とすることはできない。

18定数pivotの消去後には九本の非零Bold0残差があった。しかし81 wedge定数pivotも
消去した後の全681本について、元48成分のBold0は正確に零となることを確認する。
従ってResAT2*Bold0はこの段階で初めて全て消える。途中の九本を零と置くshortcutではない。

旧scalar operatorの全origin微分から、未知atom二次値の感度を全681本で求める。
これは二次値自体の評価ではなく、独立なsecond callbackの線形係数である。
結果にnonpure atomが残る場合、そのactual source値または別の相殺証明が必要であり、
このpreflightから数値workerを自動実行しない。

native35のrow hに出る係数はu5とcr,cl,rlだけである。quadratic core積のeven二次値では、
未知core二次値は、originがpure tauA/tauB、両nuの原点零という理由で寄与しない。
u5のeven二次値についてはdG8=N12F12のpure-tauA成分を使える。
F12の原点は零、G8の原点はtauA*tauB倍、d0=−2tauA−tauB/2なので、未知G_even2の
係数だけが単元−2となる。N2F0、d2G0の未知部分はこの成分では零。
既存二つの57 source二次cacheとphysical-F9二次adapterが全F12二次値を指定するため、
新しいgauge二次solveをせず単元−2による固定回収が可能である。dは逆にしない。

残るsource callbackの分類と全worker support上限を閉じるまで、全一次層の値やrankは主張しない。
