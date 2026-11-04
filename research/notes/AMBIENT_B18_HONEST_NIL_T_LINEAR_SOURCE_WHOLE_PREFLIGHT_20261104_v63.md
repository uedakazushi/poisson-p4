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


## Whole-cache gate の追加監査

初回感度検査では全681本に対し79種類全ての未知atom二次値への感度が零だった。
この別packetは既存source二次cacheの全葉・座標変換・C4 originの単元回収を追加で照合し、
全workerのsaved+live上限を固定する。初回感度packetとsourceは変更せず保持する。
二つの57 source解のunionは同じphysical9を共有する105葉であり、CB側を凍結済み
full9 quadratic adapterでactual CA-F9へ移してからsA=rAを課す。
追加検査もsource値そのものの評価ではない。新しいsource/gauge solveは一切許可しない。

全681列について79種類全ての旧atom二次値への感度は零だった。adapter二次項も
最終constant source組合せの全48成分で零となる。全input/cache/C4単元回収gateを閉じた。
prospective workerのsaved+live合計上限は111,057 QQ slots、60秒・2GiB・32MiBである。
新しいsource/Gauge/IFT未知数は解かず、二つの57解の保存済み二次係数だけを合成する。

最初の起動は、実在しないFULL9_FIRSTJET_REVIEWという仮のpathをmanifest段階で検出して
停止した。係数計算前であり、input-path-failure.jsonを保持する。
既存T3_ZERO_STRUCTURAL_REVIEWと実H20 source peerへ正確にbindし直した後の検査を本reportに保存した。
