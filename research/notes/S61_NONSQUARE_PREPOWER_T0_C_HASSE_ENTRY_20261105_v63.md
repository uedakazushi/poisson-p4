# 旧T0のC valuationは四：primary一回結果

2026-11-05。fresh ROOT7208 QQ/native-body/code監査後、T0専用C-Hasse workerを一回実行した。
root239.polyのa-offset4をbodyで一度だけ除いたU0=T0/a⁴を使い、runtime shiftは0。
C=100a+18e²−9e−2、a0=(2+9e−18e²)/100のHasse層を順に求めた。

**k=0,1,2,3は全係数零、最初のk=4が646項で非零**となった。
primary結果はv_C(U0)=v_C(T0)=4。逆順減算零・再現・全646係数readbackを完了した。
producerの自己一致とは別に、異なるtruncated Horner算術の独立監査を要する。

実行0.119秒、peak32832、native inner563433＋scalar/root84475 products、
入力込みRAW上界870836 bytes。34376 cells・518 bits・8MB・45秒soft/60秒hard・2GBを維持した。
guard停止やpartial prefixは発生していない。
first-nonzero-C-Hasse.txtはT0CHV001でk4、D0*100²⁰と全(e,d,integer係数)を保存する。
completed-zero-orders.txtは四つの完了ZERO層で、k≥5について主張しない。

商T0/(a⁴C⁴)や他因子・tail・GCDはまだ計算していない。
このvaluationは既知FULLU単元の因子を測ったもので、D12やT0=0の補集合を消していない。
既存の十分開・全未処理補集合・full57とquartic239の範囲を保持する。
