# 既知三一次因子によるFULLU coverのscheme飽和

2026-11-04。一回のnamed-root計算。degree4 Euclid、Gaussian、B-minor、新norm、数値点・素数は使わない。全必要unitはa方向のliteral Bézoutで証明する。

a,F,C,b1の単元性を証明し、rB=−(60t−40)/b1、rJ=(6t+4)/(5a)を作る。元のU四係数をunit*d*(d−rB)*(d−rJ)と全係数で比較する。F/Cが全体で零ならUも零であり、FULLU局所化は空。proper非単元はbase componentを保持して停止する。

H4を0,rB,rJで評価し、零の場合だけmonic synthetic divisionで因子を一つ除く。全factor identityを保存し、重複度だけ繰り返す。非零評価は単元証明を要求する。proper非単元を非零とだけ呼んで捨てない。最後にも全三評価の単元性を確認する。

この方法はKa13のirreducibilityやétale性を仮定しない。Uの全named因子はU局所化で単元だから、厳密なH=(d−r)QからHとQはその局所化で同じidealを定める。残るH_Uの各root評価が単元ならUはKa13[d]/H_Uで既に単元になり、schemeとして全FULLU coverが残る。H_U=1ならcoverは空。nilpotentや重複因子をradicalへ置換していない。

sourceは独立受理されたKa13 bridge、monicH4、FULLU係数。明示有限coverとactual Bbadへの写像・G0との組合せは親の別gateであり、このpacket単独で全域Bbad空性を主張しない。

unitcofactorの理論保存量をnative逆計算前に検査し、全体32MiB/60秒/2GBのguardを保持する。上限やproper componentに当たれば完成部分を保存して止める。

## 確定結果

承認された一回は0.9832秒で完了した。a,F,C,b1と、H4(0),H4(rB),H4(rJ)の**七つの単元証明が全て合格**した。parameter inverse、二つのnamed root、Uの四係数factorizationを全QQ係数で確認した。最後の三評価は同じ証明書を再利用できた。

除去因子数は0。残るH_Uは元H4のままdegree4で、Uはその全商代数上で既に単元となる。明示逆の有限recipeも、H_U(d)=(d-r)Q_r(d)+H_U(r)から(d-r)^-1=-Q_r/H_U(r)を三因子分掛ければ得られる。逆多項式の展開は不要である。

これは52次の明示coverをそのままFULLUへ移せるという結論で、actual badの非空・空性は結論しない。B-minor/Norm4、Gaussian、e/d-Euclid、新normは実行していない。
