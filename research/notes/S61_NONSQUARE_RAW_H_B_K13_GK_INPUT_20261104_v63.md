# Ka13上のH4/g/kの小さい共通根入力

2026-11-04。固定minor(2,3,5)は旧gそのものなので、新しいNorm4を重複計算しない。代わりに受理済みのsource g=Pi15（59項）、k=Pi20（88項）をH4で還元する入力だけを作る。d-Euclid、coefficient inverse、Gaussian、Norm4は実行しない。

sourceはPAIR_NORM_INPUT.gとADDITIONAL_MINOR_NORM_INPUT.source_kの全係数に一致する。a,e=t(a)の冪をKa13で計算し、各d係数を正確に評価する。独立な整数common-denominator raw式との差がg13で割れることを全係数で確認し、その商を保存する。

gとkはともにd次数4なので、それぞれmonicH4で一回除く。全多項式についてsource=H4*quotient+remainderをQQ[a]/g13[d]で確認する。係数先頭の未知a多項式を一つも逆にしない。

実行は0.3862秒、整数商環演算262回で完了した。gのsource次数は(a,e,d)=(4,7,4)、kは(5,8,4)。還元後はともにd次数3で、各52個の非零QQ係数を持つ。gの分子・分母の最大bit長は1536・1503、kは1638・1604。H4を含む入力全体は157個の非零QQ係数、最大1986bitである。147個のsource項の像と、二つの多項式全体のmonic除法恒等式を保存した。

次に許可を求める有限計算は、degree4/3/3のd-Euclidだけである。各係数逆元にはg13とのliteral QQ Bézoutを付け、非単元ならproper base componentを保持する。monic gcdがd-sなら、共通根coverがKa13の一個の根d=sであることをideal恒等式で示せる。その後にB21のd=s評価で全source rankを調べるのは別の予算であり、今回何も結論しない。

全norm pencil係数が零であることはH4/g/k共通根の必要coverを与えるが、全B-minorの消失を示さない。actual badの実在・空性は未結論のまま保つ。

再現は `.venv/bin/python code/s61_nonsquare_raw_H_B_K13_gk_input_20261104_v63.py --output work/GK-input-fresh`。`source-input.json` は全sourceの整数共通分母式、g13による除法商、H4による商・余りを持つ。`next-unit-Euclid-plan.json` は次段のdegree4/3/3入力と、非単元を逆にせず基底成分を保持する規則を固定する。今回の入力準備は独立レビュー済みとの主張を含まない。
