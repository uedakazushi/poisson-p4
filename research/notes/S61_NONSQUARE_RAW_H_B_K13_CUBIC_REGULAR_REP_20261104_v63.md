# 全Ka13上のcubic57に対するpositive-only regular representation

2026-11-04。入力は固定済みのactual17×24係数408個である。新しいbeta値やaの根を選ばない。既存素数32003に対し、全QQ係数分母とg13のleading coefficientが単元であることを最初に確認する。失敗した場合は別素数へ移らない。

V=Z_(32003)上でT_V=V[a]/monic(g13)は自由13である。各C_ijによる乗法を基底1,a,...,a^12で表す13×13blockを並べると、C:T_V^24→T_V^17は221×312行列になる。これは閉代数全体のrestriction of scalarsであり、aを一点に評価した行列ではない。

保存する221次minorがmod32003で非零なら、対応するQQ minorが非零で、Ka13-linear写像CはQQ上でもKa13上でも全射である。任意の剰余体へのbase changeも全射なので、全成分でrank17となる。Ka13が体・既約・reducedだという仮定は不要である。元の行列への単元scalarと既存8単元prefixはQQ上で合法であり、既受理のpointwise criterionからactual Bbad全体でfull cubic57を得る。元のGF cubic matrixやその全prefixのGF単元性を新しく主張するものではない。

nativeでは一つの221×312 work matrixでpivot列を選び、それを解放してから一つだけ221×221 minorを再構成する。partial-pivot LUをその一枚へin-placeで保存し、全48841個のP*A=L*U恒等式を元のblock係数から確認する。選択列、row permutation、Uのdiagonal、combined LU、非零detを保存する。221×221を二枚複写しない。

modular rankが不足した場合はQQ階数・full57について結論せず、次の素数・QQ巨大Gaussian・quartic計算へ自動移行しない。結果がpositiveでも、ここで証明するのはactual Bbadの全13点上のfull cubic57であり、全FULLUやquartic239・generic Poisson exhaustionへの拡張は別のsource gateを要する。

再現: `.venv/bin/python code/s61_nonsquare_raw_H_B_K13_cubic_regular_rep_20261104_v63.py --output work/K13-regular-fresh --work work/K13-regular-worker-fresh`。一回60秒・2GB・全artifact32MiB・明示150000scalarの上限。新しい係数・rank証明書の独立レビューは別途必要である。

実行は0.825秒で完了した。g13の先頭係数のmod32003像は16763で、全408分母も非零だった。regular matrixのrankは221、保存した実際の221×221 minorのdetは30542で非零。row swapは0回、全48841個のPLU係数恒等式とdetのpivot積を検証した。nativeの最大scalar slot上界は74269、全artifactは約311KBだった。

従って保存係数と既受理のsource gateの下で、**actual Bbadの全13幾何点でfull cubic rank57が成立する**。Bのrankは落ちるが、他のpure/genuine方向がfull cubicを救っている。このpacketではwhole FULLUの残る別stratumやquartic239を調べず、追加計算も自動実行しない。
