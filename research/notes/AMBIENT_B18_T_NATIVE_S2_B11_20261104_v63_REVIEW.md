# native S2 B11：独立source・有効代表・消去条件監査

2026-11-04。個々のactual native X2を求めたという主張をしない。
出力は零originのsource結合に使える `(L0+L1) S_actual_B11` である。
作者producer/helperは呼ばず、既存独立tensor engineと別の疎B11/C44積で検査する。

実B11は基底 `1,r,r²,z0,z1,z4,z5,z0r,z1r,z4r,z5r`、idealは
`(r³,zi zj,zi r²,z2,z3,z6)`。実H20の全400積をこのmonomial積と照合し、
全105保存rootを正確に写す。実pair `IA=(z7,r³)`, `IB=(z7+c r²,r³)`の全積が消えるため
環写像は全次数で正当であり、同じ長さや二次jetだけから推定しない。
両monicを実H20値から運び、C44の全1936基底積を確認する。

全source DAGをB11で評価し、保存物理150 trace pairに一致することを検査する。
固定30系の150補正を実constraintへ代入し、全30 scalar B方程式零を確認する。
補正はm_B²に属しm_Bが殺すので、既受理origin J30による線形式とその一意性が使える。
補助etaのscalar係数を確定してから同じ係数をphysical tauAへ埋める。
全actual14・κ9・G8をC44で評価し、native Sの全2109係数と比較する。
新しい105求解・物理source jet・N2/R12/T2の計算はない。

scalar Lの0/1係数は、まず既受理C36で `Res4(N^t)` を取り、そのscalar係数だけをB11へ運ぶ。
native Nをtruncateしてsecond monicでRes4する方法は、特定のL2を誤って残すため採らない。
全128列の有効代表6876係数を照合する。

source-specific L2消去は以下の通り。old640のorigin像rank4、native128のrank23、合計rank27なので
二像は直和である。全741と620の結合を全56成分で再評価し、old/nativeが別々に零と確認する。
L0が単位だから `S0 c_native=0`、従って積の未知項 `L2 S0 c_native` は射影前に消える。
物理sourceに入るS2やmonic二次を零とはしない。

m1消去について、620結合の全56成分で `x1=Xp0 a1` と既知Pi3の全464二次関係消去を再確認する。
`z=x-Xp a1`は二次から始まり、`(mz)2=m0 C2`。未知a2はm0とpiがXp0を殺すため消える。
EulerはR線形でないが、そのcommutatorはKoszul cycleで、一次は
`H2 K(z1,...,z7,0)=H2 K(z7)⊕H1 K(z7)=0`によりliteral d3一次の境界である。
五targetへ運んだ差は既知d3二次spanに入りPi3が殺す。これは三class商だけの等式であり、
未射影のSchur二次そのものとのliteral等号ではない。

個々のX2、全620classの値、T長さ、K torsion、元深さは未決。
専用checker `code/ambient_B18_T_native_S2_B11_review_20261104_v63.py` と同名REVIEWに保存する。
