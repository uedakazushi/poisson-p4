# cubic57：projected15 filtration と二行heavy係数の有限preflight

2026-11-04。元57×91 source の normalized43,49（actual135 cubic20,61）の
H-linear七係数だけを新たに抽出した。追加core点、rank計算、GB、四次係数展開は使わない。
個々の HZ/HW係数は既知 a,d,F,B,C,J の積の単元ではなく、単独係数のunit gateは不成立。
両係数は非零であり、これを同時零と扱わない。pairの明示的組合せは別のUNIT_PAIR packetで検証する。

FULLUではsix-quad degree3 product imageはrank41 subbundle、ZW制限のfour-quad product imageは
rank19 subbundleである。quadにHは無いので、degree3商Bは pure43⊕H*light7 のrank50 bundle。
ZW制限の商Vはpure16⊕H*five のrank21 bundle。
従って0→T29→B50→V21→0はexact。degree3のこの議論にはU78、CI、ell regularityを使わない。

compact source identityによりdelta*(G61−G20)はZW制限後にquad productsへ入る。
ここでdelta=−64000C²はFULLU unit。
既存のrestricted36×40の保存nonzero34minor Delta34が非零の開では、projected genuine16は
rank15で、そのkernelはprimitiveなG61−G20 lineとなる。
従ってこの開でfullgenuine16がquad product41から独立であることは、
delta*(G61−G20)のB50内のclassが各fiberで非零であることと同値。
quad productsはHを持たないため、HZ/HW coefficientへのunit functionalがあればこの非零性を示せる。

実sourceの七係数差は、oldlight=(x0,x1,x2,x3,x4,Z,W)へ戻すと
最初の五係数がすべて零で、kZ*Z+kW*Wの形になることを直接検証した。
個々の既知unit testはfalseであり、このpreflight単独では全域の非零classを証明しない。
数学データはcertificate.json、実行時receiptはprepublication-receipt.jsonとして保持する。
全FULLU cubic57、fullquartic239、knownkernel38/missing42はこのpacketで変更しない。

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_cubic_projected15_heavy_preflight_20261104_v63.py --output work/marked-square-v63/cubic-projected15-heavy-REPLAY
```
