# 孤立無限遠axisの実polarization idempotentとcharted sheaf分離

2026-11-03。直前のcolon recipeに、元134式の実際の45変数quadratic chartから作るidempotentを加える。未知アフィンの全schemeを保持し、証明済みの孤立無限遠因子だけを分離する。全w-torsionを消したとは主張しない。

## 1. 実sourceと全midpoint identity

保存axis pointで非零のq一座標をそのsample値へ固定するprojective chartを取る。他の44qとwが45座標xであり、axis point pはp_w=0を持つ。元134式を全てq次数2へwで斉次化して、このchartへ直接代入する。producerは全134のF_i(p)=0と、全QQ係数の

    F_i(x)-F_i(p)=J_i((x+p)/2)*(x-p)

を検査する。保存45行は全て元134式に属し、そのd=det J_selected(p)は

    -124102513209883356037120000000000000000/59049 != 0

である。全134 chart式と45×45 linear midpoint matrixの係数を保存する。

## 2. 実ideal分解

D_p(x)=det J_selected((x+p)/2)、e_p=D_p/dとする。selected45 identityにadjugateを掛けるとD_p(x)(x-p)∈I^45を得る。従って

    e_p*(x_j-p_j)=0 inA_patch

である。e_p(p)=1なので、多項式の評価差を(x-p)で因数化すればe_p(e_p−1)=0となる。さらにe_p f=f(p)e_pであり、e_p因子は実際にKという一つの被約点である。従って

    I=(I+(1-e_p)) intersect (I+(e_p)),
    (I+(1-e_p))+(I+(e_p))=O

がスキームとして成立する。p_w=0から**w*e_p∈I**である。D(w)ではe_p=0なので、補因子I+(e_p)は元Iと正確に同じになる。未知の単純アフィン点も、正次元や非被約な未知構造も削らない。determinantの次数は≤45で、巨大展開はしていない。

## 3. generic finite etale axisへのtrace

先のsheaf-colon recipeは、global homogeneous45Jacobian minor Deltaを使って

    I_axis=(I+(w)):Delta^infinity

というfinite etaleかつclopenなgeneric axis factorを与える。Euler関係によりDelta≠0なら、minorで除いたq座標も非零となる。従ってこのfactor全体を同じprojective chartで扱える。betaを必要な非空開へ縮め、有限etale代数をE_axisとする。

E_axis上の普遍点pへ上の実quadratic polarizationを適用する。d(p)はE_axisの単元となり、そのnormによる追加局所化はbetaだけである。separable closure上で異なる点のe_pは、どれかの座標差により互いに直交する。従って

    e_axis=Tr_(E_axis/K)(D_p(x)/d(p))

がgeneric axis factor全体のidempotentとなる。全e_pでp_w=0なのでw*e_axis=0であり、補ideal I+(e_axis)もD(w)上で元Iと一致する。これは有限etale代数の存在と明示的なdeterminant/traceによる有限recipeである。E_axisの次数、basis、traceの展開係数は計算していない。

chart外にはaxis factorがないのでidempotentを零で延ばし、projective Yのsheaf idempotentとしてglueできる。これは `I_rest=I:I_axis` による補sheafと一致し、全未知アフィンを保つactual ideal coverとなる。generic Kと共通非空beta開での主張で、全special fiberでの一様な分解は主張しない。

## 保存と範囲

`code/s61_affine_axis_polarization_split_20261103_v63.py`、`reports/S61_AFFINE_AXIS_POLARIZATION_SPLIT_20261103_v63/`。実QQ sampleで全134 midpoint coefficient identitiesとdを再生し、generic axis traceのproofをここに与える。全w飽和、残りraw成分のpurity、未知アフィンの空性は未完了である。
