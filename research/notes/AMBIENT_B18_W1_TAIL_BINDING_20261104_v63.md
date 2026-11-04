# Old49からactual W1 tailへのliteral単位Hensel表示

2026-11-04。旧218式の全交叉を、現在のphysical G19切断とw=1へ実際に戻した。
別々のnormal formを同一視せず、旧49座標から元physical150係数への全mapを保存した。
Jinfのsyzygy、T3のpresentation、元B18原点の全D depth/nilは未計算である。

## 1. 旧49座標の同じambient embedding

旧変数はz0..z34,bb,cc,v0..v9,ss,aaで、中心はv4=v7=1、他は0。
新しい形式変数はこの中心からの差として取る。
旧critical beta sliceはss=s,aa=p²,bb=p³,cc=spであり、保存された全beta24多項式を
これら四変数へ戻した。固定adapterの有理左逆でbeta16を得る。
q45は旧selected10をv、remaining35をzとしてそのまま戻す。
これでS61の61座標が全て指定される。

元S61 certificateの55 eliminated係数を使って150係数を再構成した。
ここはEuler-RREF representativeなので、そのままphysical cutを読まない。
divergence-free representative

```
P_can = P - (X*div(P)^T-div(P)*X^T)/5
```

へ正確な線形mapで戻してから、旧W1固定GLの逆を使う。
原点の全150係数が既保存のE(1,0)と一致することを確認した。

## 2. 実際の17次元group sectionと全20切断

新W1 frameで六value方向E22,E44,E32,E42,E31,E41と四perspective方向E0jを選ぶ。
各方向にはI+tEijを使い、diagonalの場合の逆元だけ1/(1+t)を使う。
旧Sp4 orbitのselected7方向には

```
C_A(t)=(I-tA/2)^(-1)*(I+tA/2)
```

という有理symplectic Cayley sectionを使う。
全逆元の分母は原点で非零である。
17因子の積を固定GLで元frameへ運び、再構成した150係数へ作用させた。
このpure rational DAGには不明なbranchのpullbackや非単位除算は無い。

common G19の19 physical cut rowsを基点値へ戻し、physical row140を1とする全20条件をcとする。
前packetのselected17条件に関するgroup17 Jacobianは、今回のliteral DAGのADでも
正確にdeterminant **32**。従って17未知数は旧49変数の一意なpointed algebraic power series。
このunit Hensel系を全保存し、解のjetは計算しない。
残り三条件はpositions14,18,19であり、最後がw=1である。
それらのliteral nodeも保存した。

従って候補を省略せず保持したtail表示は

```
QQ[[旧49中心差,group17]]/(旧218式,全20 physical切断)
```

であり、unit17未知数を除けば旧49の218式に三つの明示的なalgebraic関数を加えた表示になる。
physical R0の九座標も、同じ150出力のrows
27,108,109,112,113,117,118,119,61から読み戻した。

## 3. 型の一致ではなくactual germを与える理由

旧W1認識は、固定GLで移した実テンソルが旧critical sourceの
p=s=r=0,k0=1,k5=1/3に文字どおり一致すること、両成分のfiber blockがdet -72で
可逆であること、および旧49表示が各実image kernelの全kernelであることを保証する。
value/perspective10方向とcritical beta orbit7方向の同時形式積も既受理である。
今回の17 sectionはまさにこれらの方向を接線とし、各blockが単位なので、同じ形式積の
別の有理sectionを与える。common base-imageの五つのbeta graph関係は両成分idealに入り、
除いた後の49 ambient embeddingも上記の同じmapである。

group作用はCA、CB、Nとその実idealを保つ。
したがって二つのactual kernelの和も保たれ、旧218式との対応は非被約な全Dに対して成立する。
接空間rankや個別Ac型だけから交叉を推測する議論ではない。

現在のw次数表示についても、元sourceのzero48系とpositive5系は実normalizer chartであり、
正次数の出力は正確な有限w多項式である。W1では全normalization fiberが一つの幾何点である
という既受理のproper-fiber結果と、同じsource chartのunit条件があるため、w=1の像が
別の小さなsource枝だけを拾う問題は無い。
従ってこのliteral切断環は、現在の安定module H_A/Jinfのactual ring表示を与える。

## 4. 三つの残余条件のregularityとtailの次元

元のG19 section（tensor scalarとGL18方向）についても、同じE1点でcommon cut19への
Jacobianを独立に作った。そのdeterminantは **-32**。
ambient canonical係数空間でgroup作用とphysical切断のunit IFTが使えるので、
D、N、そのnilmoduleに対して同時のG19形式積が成立する。
したがって19切断はこれら全てにregularである。

切断は既存gradingのtorusに保たれ、wはweight1、w(E1)=1。
wを逆にすると正確に

```
D_perp[w^(-1)] = (H_A/Jinf)[w,w^(-1)]
```

となる。同じことがNとnilmoduleにも成立し、w-1はさらにregularである。
この20個のregular切断を、別のG17形式積の17単位未知数を先に消す順へ並べ替えると、
残る三条件が旧49のD、N、nilmoduleにregularであることが従う。
local有限moduleの最大ideal内のregular sequenceの順序交換を使う。

旧49のDとNはdimension11、depth4、nilmoduleはregular4次元上free1。
従って今回のtailはDもNも **dimension8、depth1**、nil supportは1次元となる。
被約Nの8次元成分を消してwhole ringを一変数の二重環へ置き換えない。
これらはw=1 tailの結果であり、元w=0 B18原点のdepth19/20ではない。

## 5. Nilを同じR0上で指定する

旧受理済みのnil generatorとannihilatorをそのまま保持した。

```
F_old = 2*aa+v0+z25,     F_old²=0,
ann(F_old)=旧45本の線形式。
```

annの自由変数は(v1,v4,v5,v7)。group17未知数をunit系で除いた後、
残余三条件のこの四変数上の微分は

```
[[0, 0, 1, 0],
 [0,-1, 0, 0],
 [0,-1, 0,5/2]].
```

rank3 minorは **5/2**。従ってnil support curveは実際にregularであり、
その接線はv1方向である。
この接線に沿ったphysical R0九座標の速度は、row117だけ-1、他は0。
したがってphysical R0 row117がcurveの正確なuniformizerとなる。

さらに、旧ann45を課した全20切断と
`physical_row117 - T=0`からなる21未知数のunit入力も保存した。
未知数はgroup17と旧free4で、Jacobian determinantは **-80**。
これにより全九physical R0関数をTのalgebraic power seriesとして一意に指定できる。
R0からcurveへの写像は完備局所環として全射であるから、nil tailは

```
(F_old) ≅ (R0/P_nil)*F_old,    R0/P_nil ≅ QQ[[T]],
```

という同じR0上のcyclic moduleである。
P_nilはrow117以外の八座標を保存した21-unit inputのreadback関数へ等置するregular ideal。
curveに制限するとF_old自身の関数値は0になるが、moduleのfree basisとしてのclass F_oldを
0と同一視していない。

三条件がNにもregularなので、旧nilの完全列は切断後も正確であり、nilmoduleは失われない。
tailの被約商は旧218式にF_oldを加えた同じ切断環である。
この表示は同じR0上のexact sequence

```
0 -> R0/P_nil * F_old -> H_A/Jinf -> N_tail -> 0
```

とそのactual ring extensionを指定する。
F_oldのH_A=(1,tau_A)二行lift、N_tailのR0上の最小module presentation、Jinfの全syzygyは
まだ構成していない。従ってT3のrelation moduleが完成したとは主張しない。

## 保存・実測

`code/ambient_B18_W1_tail_binding_20261104_v63.py`、
`reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63/literal-tail-binding-DAG.json`。
全10516 scalar node、旧49→S61→canonical150、groupとinverse、218式、20切断、
17-unit系・残余3条件、R0九関数、旧nil/ann、nilcurve21-unit系を保存した。
一回0.686322秒、最大RSS121664KiB、証明書169982bytes。
Newton jet、GB、新点、source再生成、opaqueなbranch pullback、p逆元は無い。
