# Ka13全体でのH4/g/k共通根の互除法

2026-11-04。入力は固定済みのK13_GK_INPUT（report SHA 632e472246977ecf915eb6eabb42c8949caee98f05f1ab2ba2f5c514920e3282）のH4・gR・kR、d次数4/3/3である。Ka13=QQ[a]/g13を体とも既約とも仮定しない。先頭係数を逆にするときは、primitive整数分子pについてS*p+T*g13=R≠0を保存・全係数確認し、その後に逆元を使う。真の非単元なら全基底成分を保存して止める。

各単元証明を先にcheckpointへ書き、その後にmonic化・d除法・source cofactorの計算を進める。除法A=Q*B+Rと単元によるmonic化は全係数で検証する。この連鎖は生成idealの両方向の包含を保存する。最後にgcdが各入力を割ること、gcdが三入力の明示Bézout結合であることを検証する。全記録は`checkpoint.json`、最終cofactorと根は`result.json`に保存する。

gcdがd-sなら(H4,gR,kR)=(d-s)がKa13[d]上のschemeとして成り立つ。これは全B-minorの消失ではない。B21への代入、別minor、Norm、数値点・素数・GBは許可範囲に含めず、このゲート後に自動実行しない。

入力は157個のQQ係数・最大1986bit。実行上限は60秒・2GB・32MiB・明示多項式150000項、先頭係数単元照会は高々9回とする。逆元証明のHadamard上界と保存済みartifact量を各照会前に確認する。上限・非単元・ゼロ先頭係数の場合は完了と呼ばず、checkpointを保存する。

再現: `.venv/bin/python code/s61_nonsquare_raw_H_B_K13_d_Euclid_20261104_v63.py --output work/K13-d-Euclid-fresh --work work/K13-d-unit-fresh`。新しい入力点を選ぶ操作はない。この計算の独立レビューは別途必要である。

結果は1.309秒で完了し、gcdはKa13全体でd-sとなった。先頭係数の単元証明は4件、全多項式nodeは15、商環整数演算は567回、最大係数bit長は3458だった。sは13項、分子最大1531bit・分母1490bitである。明示source Bézout cofactorのd次数はH4に対して1、gRに対して2、kRに対して零多項式であり、係数はそれぞれ26・39・0項。三入力すべてのgcd除法の余りが零であることを全係数で検証した。artifact総量は約1.67MBである。

元のsource g・kについても、固定入力のg=H4*qg+gR、k=H4*qk+kRによりidealは同じである。保存したd-s=A*H4+B*gR+C*kRから、d-s=(A-B*qg-C*qk)*H4+B*g+C*kという明示cofactor recipeを得る。これは新しい逆元を使わない。ここでC=0である。したがってH4/g/kの全共通根schemeはKa13であるが、残りのB条件の成否はまだ別問題である。
