# 固定minor(2,3,5)は既存gでありNorm4は既に零

2026-11-04。承認された固定minorのsourceを先に束縛した結果、重複する既知零normだと判明した。新しいNorm4 arithmeticや別minorは実行しない。

lexicographic B35でrow(2,3,5)はindex26。そのlocalized generator mapは15で、既存PAIR_NORM_INPUTのgそのものである。prepared minorの全59係数が、このgの一つの非零QQ倍と厳密に一致する。original source minorへのrow/column preconditioningは既知FULLU単元だけで、現在のH4 coverではU自体が単元と証明済み。

旧source normはknown-unit factorsを除いたoldhを持ち、oldh=G*h1を全QQ係数で確認する。現在のKa13上ではh1(a,t)=0。したがってNorm(g)=0、固定prepared minorのNormもQQ単元の四乗倍なので零である。これはGをglobalに逆にした推論ではない。h1=0からoldh=0へ掛け戻すだけである。

このsource gateは、候補を支持する新しい実験ではなく、既に分かる零を再計算しないための記録である。actual badの存在・空性は未結論。新しいminorへ自動切替せず、親へsource lineageを報告した。
