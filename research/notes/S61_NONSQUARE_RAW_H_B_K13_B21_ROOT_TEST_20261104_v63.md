# Ka13共通根上の全B21と有限minorゲート

2026-11-04。受理済みのB21係数artifactと、新しく固定したd-Euclidの根d=sを使う。Ka13=QQ[a]/g13の全体を保持し、体や既約性を仮定しない。7×3の全21成分はd次数3以下なので、Horner法の乗法は63回である。先に全値とB_ij(d)=(d-s)*Q_ij(d)+B_ij(s)の係数記録を保存する。

その後、行0..6の辞書式3部分集合で高々35個の3×3 minorを検査する。最初の非零minorについてだけg13との整数Bézout証明を付ける。単元なら全Ka13上のrank3を証明する。真の非単元ならその基底成分を保存して止める。35個すべてが零の場合、このpacketは候補を保持するだけで実在を主張しない。

source bindingは次の通りである。MINORS_COSTの全21恒等式はrowunit_i*Bprepared_ij*colunit_j=Boriginal_ij*D_jで、row/columnの因子はQQ定数とFULLU=adFBCJの既知単元のみである。したがってFULLU上の最大minor idealは不変である。B21入力の全84係数像はK13_SOURCE_BRIDGE_REVIEWで独立受理されている。GのR/Iでの単元性、norm pencilのsource membership、Ka13のe=t同一視、今回の(H4,g,k)=(d-s)により、Ka13からactual R/Iへの全射がある。従ってこの候補上の単元minorはactual bad schemeの空性を与える。Gが原環の単元だと仮定する議論ではない。

この計算は元の22条件・source openの全検証が済んだとの独立主張をしない。全minor零の場合のactual存在判定は、その全source gateの追加確認を待つ。単元minorの場合でも、full cubic57・quartic239・generic Poisson exhaustionは自動的に従わない。

上限は一回60秒・2GB・32MiB・明示150000項、単元照会は高々一回である。新しい点、素数、GB、Norm、別の巨大行列は計算しない。再現: `.venv/bin/python code/s61_nonsquare_raw_H_B_K13_B21_root_test_20261104_v63.py --output work/K13-B21-root-fresh --work work/K13-B21-unit-fresh`。新しい結果の独立レビューは別途必要である。

実行は0.627秒で完了した。全21値はそれぞれ非零の13係数で、合計273個のQQ係数、最大分子1672bit・分母1606bitである。その後の全35個の3×3 minorはKa13で厳密に零となった。単元照会は0回で、基底成分は一つも捨てていない。商環整数演算679回、2×2 cofactorを含む最大係数bit長1780、artifactは約2.83MBである。

したがってこの有限候補を空とする試みは成功していない。保存した結果はrank≤2を示すが、このpacket単独ではactual Bbadの実在を主張しない。次のsource gateは、元の22生成元の全体とFULLU開条件を同じKa13の(a,t,s)で確認し、既知のcover全射と両方向を照合することである。これはここでは未実行であり、別の計算・レビューを要する。full cubic57の真偽もこの結果だけでは判定しない。
