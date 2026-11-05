# T1の既知C四乗商4811項と全multiply-back零

2026-11-05。ROOT_SOURCE/code/budget独立監査後、一回だけ実行した。
Q=T1/(a³C⁴)を有限Laurent polynomial-part operatorで計算し、
**全係数でC⁴ Q_integer=100²⁰ (T1/a³)_integerを確認した**。
負のa冪remainderを仮定で捨てず、全multiply-backにより零とした。
各slabも途中の高次e項を全て保持し、完成後にだけ想定支持外の零を確認した。

Qは4811項。7595項の入力に比べ項数が減った。
sourceのwhole81776 cells・全整数518 bits・RAW10.61MB・45秒soft/60秒hard・2GBを維持し、
実行0.121秒、peak56112、polynomial products467824、scaled/root products83013、
入力込みRAW上界1454299 bytesだった。guard停止やpartial prefixは生じていない。

出力はquotient-integer-common-D.txtのT1C4QV001で、global Dq=D1*100²⁰と全packed key／整数係数を保存する。
全4811係数をreadbackした。factorizationやscalar content normalizationはしていない。
全multiply-backはprimaryの厳密検算であり、別long-division等による独立監査は別gate。

v_C(T1)=4の独立Horner peerを合わせれば、QはCでこれ以上割れない。
除いたa³C⁴はFULLU既知単元なのでFULLUでの零集合は変わらない。
この計算はGCD、他の因子、他tail、strict cover enlargement、不良集合の除去や
全FULLU cubic57を主張しない。次のSOURCE予算も別gateである。
