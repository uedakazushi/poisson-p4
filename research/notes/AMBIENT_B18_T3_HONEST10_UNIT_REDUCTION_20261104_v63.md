# Honest10 presentationの単位Schur消去：残る最小generatorは9個

2026-11-04。既構成の正確な10行presentationを原点で評価し、見つかった単位関係を
有限DAGで消去した。T3の零・非零、colon、全Dのdepth/nilはまだ未決定。

## 正確な原点係数

凍結済みの全186 DAGに対する四つのg商原点値0のoverrideを保持し、
honest10構成で追加した792 coordinate商には保存済みの正確な原点値を使った。
新しい全scalar nodeを有理数で評価し、全inverseの分母が単位であることを確認した。

36Koszul+2034 Q[3] liftsの10×2070原点行列のrankは **1**。
top12を追加したD4行列、および36Koszul+2046 Jinf liftsの原点rankも1で一致した。
top12自体の原点行列は零。
元の二行matrixでは見えなかった単位関係が、一つだけ現れた。

## 一つの厳密な単位消去

使用したrelationは0-basedで全2070列中1389列、すなわちQ[3]1353列である。
これはweight1のactual physical graph row52を係数1で用いた列。
original free generator5、すなわちf_CA60*tau_Aの係数は原点で **10/3** なので単位である。

current relation vectorをc、消すrowをi、単位u=c_iとする。
任意のrelationやtop vector vを残る各row kで

```
v'_k = v_k - c_k*u^(-1)*v_i
```

へ移す。この式はrelation cによってgenerator iを解いた正確な商moduleの同型であり、
有限jetによる消去ではない。uの逆元は非零原点10/3を持つunit inverseだけを使う。
全relationとtop12を同じ射影で移し、元の10-generatorから残る9-generatorへの
明示的な9×10射影も保存した。残るoriginal generator上で射影が恒等となることも確認した。

その後、全残余relation係数の原点値は零となり、次の表示を得た。

```
D3 = coker(R0^2069 -> R0^9),
D4 = coker(R0^2081 -> R0^9),
2081列 = D3の2069列 + 同時に移したtop12列。
```

明示DAGとして恒等零になった列以外は捨てない方針で検査したところ、
この段階ではそのような列は一つも無かった。
この2069/2081は完全な関係の本数であり、最小relation数とは主張しない。
全relationが最大idealに入るので、**mu_R0(D3)=mu_R0(D4)=9** は正確に従う。

T3はD3内の変換済みtop12の像である。
top12が原点で零であることはT3=0を意味せず、T3のpresentationを得たことにもならない。
次に必要なのはこの2069関係に対するtop12のrelation moduleまたは実membership証明である。
期待した二、三generatorへの単位消去は起こらなかった。

## 範囲と保存

古いnonunit r/g関数は同じactual physical baseのscalar coefficientとして保ち、
restriction・再除算・未知微分の零置換をしない。全Koszul関係も含めた元presentationから
正確な単位同値変形のみを行ったので、ambient二行matrixのTor問題を無視していない。
GB、jet、新点、source再生成、非単位逆元は使わない。

Source `code/ambient_B18_T3_honest10_unit_reduction_20261104_v63.py`。
`reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63/unit-reduced-presentation-DAG.json`に
全498311 scalar node、792原点値、全原点行列、unit pivot、射影、全残余matrixとtop12を保存。
証明書9964518bytes。一回1.298746秒、最大RSS266044KiB。
