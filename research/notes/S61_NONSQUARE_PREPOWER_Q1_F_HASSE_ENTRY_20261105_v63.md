# 正規化C四乗商のF valuationは零

2026-11-05。fresh ROOT source/code/budget監査後、許可されたF-Hasse workerを一回実行した。
Q_C=T1/(a³C⁴)、F=500a+90e²+3e−24。
ROOT独立longdivision出力の4811整数係数・D197をbyte同一で使い、
R_F=24−3e−90e²、a0=R_F/500へHasse代入した。

**最初のk=0が665項で非零**となったので、その場で停止した。
primary結果はv_F(Q_C)=0であり、Fでこれ以上除く商は作らない。
逆順減算で零へ戻して再現し、全665係数をreadbackした。
独立の異なる算術Horner検証は別gateであり、三passの自己一致と同一視しない。

実行0.0345秒、peak owned27530、native inner products190503、scalar/root products29651、
入力込みRAW上界608643 bytes。
固定28708 cells・hard518 bits・8MB・45秒soft/60秒hard・2GBを維持した。
guard停止やpartial prefixは発生しなかった。

出力first-nonzero-F-Hasse.txtはQ1FHV001、k=0、一個のD0*500¹⁶と全(e,d,integer係数)である。
completed-zero-orders.txtは空で、未評価のk≥1を零と主張しない。
この結果は既知単元Fの因子が無いことだけを示し、T=0補集合、不良集合、
他因子・tail、GCD、全FULLU57や四次239には新しい主張を加えない。
