# Q₁でのfirst CB branchとCA：actual closed transverse CI4

2026-11-03。actual D全環ではなく、Q_bのb=1での**一つのsmooth CB branchとsmooth CAの実intersection**を調べる。
先行の非対称な接rank36から共通接空間は22次元であり、実GL×b paired-image S'はsmooth20である。

## 実58 source座標と36 implicit方向

CBは実19fiber+ratio+flag9、CAは実20fiber+flag9で、それぞれsmooth29 sourceである。
同じ物理点へ戻した全150係数differenceのJacobianはrank36。
共通source-kernel22のうち、実GL×b接方向をmatched pairとしてsourceへ戻し、20個のS'方向を取る。
残り二方向を選び、36個のpivot source座標と合わせた全58×58 frameのdet非零を確認する。
36×36 target/source Jacobian blockの正確な逆行列も保存した。

新しい座標をbase p20、implicit y36、normal(r,z)とする。
pは実S'上でétale座標となる。36式のformal inverse function theoremでyを**全次数で**消去できる。
閉fiber p=0の全residualは二normalの形式級数であり、一次項は零である。

## 物理的な全二次係数

全150行のsourceを実5×5指数行列とその逆座標で代入して、二normalの全二次係数をpolarizationにより計算した。
CBでは常にg(t)=mu(t)f(t)を保持し、t² delta_mu delta_fを落とさない。
今回選ばれた二normalのratio speedは0であることも明示する。
36次の単元blockでimplicit二次補正を施すと、全150×3 obstruction行列のrankは2。
独立なrow63,68（0起算）の二式は厳密に

    Q1=6(r−z)²,
    Q2=−23r²/15+18rz/5−9z²/5。

Q2+3Q1/10=4r²/15であり、したがってそのidealは(r²,(r−z)²)。
二式の一次multipleが全三次四monomialを張る4×4 determinantは64/25である。
これは有限fieldのrankや少数方向sampleではなく、全Sym²(two-normal)のQQ係数を持つactual obstructionである。

## なぜ全closed fiber idealが二quadricsに等しいか

二つのactual residualをF1=Q1+O(m³)、F2=Q2+O(m³)とし、m=(r,z)とする。
三次multipleの単元性からm³⊂(F1,F2)+m⁴。
有限module m³/(m³∩(F1,F2))へのNakayamaによりm³⊂(F1,F2)。
従ってQi=Fi−O(m³)も(F1,F2)に入り、逆包含もm³⊂(Q1,Q2)から従う。
全150 residualの残りは、その二次項をQ1,Q2のconstant combinationで除くとO(m³)になるので、既にこの二式のidealに入る。
従って**実closed transverse fiberの全completed ideal**は正確に

    (r²,(r−z)²)

である。単に2jet idealをactual ringと呼び替えたものではない。
closed fiberはlength4、HF(1,2,1)、m³=0のGorenstein complete intersection。
座標(r,r−z)でC[[u,v]]/(u²,v²)となる。

complete Nakayamaから、実completed CB1∩CAはregular20 base上にfiniteで高々4generatorを持つ。
この一点のclosed lengthからgeneric length4、相対flatness、全relative ringを推論しない。
actual Dのnil module・primary decomposition・第二embedded primeも未決である。
generic QQ(b)のquadratic rankと有限rank下限を別に検証することが次の課題である。

再現:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CB1_CA_transverse_CI4_20261103_v63.py --output <new-output-directory>

全58frame、36inverse、二normal、全150raw/corrected quadraticを保存する。
GBなし、120秒・2.2GB。
