# Ka13の導関数単元証明

2026-11-04。actual B-rankbad algebra=QQ[a]/g13、長さ13は独立受理済みである。本packetはそのsource bytesを変更せず、g13と導関数g13'について一回だけextended gcdを計算する。既約因子探索やKa13を体とする仮定は不要である。

導関数の整数contentだけを除いたp=g13'/cについて、保存する整数恒等式S*p+T*g13=R≠0を全係数で確認する。さらにS*g13'+c*T*g13=c*Rへ戻す。gcd=1ならQQ[a]/g13は有限étaleであり、代数閉包上で相異なる13個のreduced pointに分裂する。従って独立受理済みのactual scheme同型により、actual Bbadも幾何的にreducedな13点である。

これは各点がQQ有理点という主張でも、g13のQQ既約性でもない。指定したBブロックのrankbadについての結果であり、full Poisson cubic matrixのrank低下やPoisson schemeのnilpotent structureへは拡張しない。

再現: `.venv/bin/python code/s61_nonsquare_raw_H_B_K13_squarefree_20261104_v63.py --output work/K13-squarefree-fresh --work work/K13-squarefree-worker-fresh`。60秒・2GB・32MiBの一回計算で、新点・素数・factorization・GBを使わない。

結果は0.494秒でgcd(g13,g13')=1。導関数・g13へのcofactorは13項・12項、最大3692bit、非零整数定数は3779bitである。primitive導関数への恒等式と、元の導関数へ戻した恒等式をともに全係数で検証し、約66KBに保存した。新しい証明書の独立レビューは別途必要である。
