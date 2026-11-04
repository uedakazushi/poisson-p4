# 五点・射影S5商・stable level 2：独立監査

2026-11-04。対象は frozen `Q_FIRSTJET_FIVEPOINT_S5_REDUCTION_20261104_v63`
と、その新しい `Q_FIRSTJET_FIVEPOINT_VECTOR_DESCENT_20261104_v63` companion。
元W16のSp4不変**関数体**Fを扱う。有限群不変式探索、β点、高次数probe、GB、
作者producerの再実行はない。一回の小さい整数行列検査と以下の幾何学的証明を用いる。

## 五点と5次元の線形表現

一般U(2,1),V(1,2)はQ=P1×P1上で横断的に交わり、交点数は5。
Koszul列をO(2,1)でtwistすると、両端の必要なH0/H1がKünnethで零となり、
h0(I_Z(2,1))=1。同様に(1,2)も1。従って交点集合から両projective形式が
一意に回復できる。sourceと一般unordered五点configurationは既約10次元で、
char0で一般fiber一つなので双有理である。ruling交換とU,V交換は同じ作用であり、
全Aut(Q)を割るときにlabelやorientationの余分な二重被覆を残さない。

P3で五つの一般ordered pointsはprojective frameである。smooth quadricsと
frameのincidenceをPGL4で割ると、固定frameを通るquadric線形系P(V5)になる。
逆にquadricを固定すればconfiguration/Aut(Q)。labelを外す作用はS5である。
frame vectors f0,...,f3=e0,...,e3, f4=−Σei は、S5の**真の線形Std4作用**を
持つ。したがって

    B=C(projective pair)^Aut(Q)=C(P(V5))^S5.

V5はSym²(Std4*)→permutation5の核。直接の係数計算では、四つのdiagonal
係数が零、六off-diagonal係数の和が零となり、次元5である。
独立checkerは作者のsymbolic substitutionを使わず、整数linear-form積から
全120個の5×5作用を再構成する。保存四生成行列、Coxeter関係、全120×4の積、
class characterとprojective faithfulnessを確認する。characterは
`(5,1,1,−1,1,−1,0)`、内積1でV5=[3,2]。非恒等元はprojectively scalarで
ないのでその固定locusは真の閉集合であり、有限和の補集合でP(V5)のS5作用は自由。
これにより一般configurationのAut(Q) stabilizerも自明となる。

## 二倍率の降下：Brauer障害を残さない証明

既受理Spin4切片・追加normalizerの完全な同定を用い、非零affine pairの空間を
`X=L_U^××_Y L_V^×`, `Y=P(U)×P(V)` と書く。
中心Z=μ2²はU,Vそれぞれの独立な符号なので、**実際のsquare map**により

    X/Z = {両成分が非零な L_U²⊕L_V² の開部分}

となる。二つのsquared tautological lineの直和は、中心が自明に作用するため
PGL2×PGL2でhonestly linearizedである。追加成分τはsummandを交換し、その
二乗は元の中心なのでsquare上で恒等。従って全Gext/Z=Aut(Q)に対する
**真のrank2 vector-bundle linearization**を持つ。

一般自由なAut(Q)作用のRosenlicht quotientの一般点では、Yはtorsorである。
このvector bundleをtorsorに沿って降下すると、B上の普通のrank2 vector space
を得る。両成分非零の開はdenseで関数体を変えない。従って

    F=C(X)^Gext=B(s1,s2).

rulingを個別にB上のP1としてsplitする仮定は不要である。交換されるlineの
**直和そのもの**が降下する。spin liftの2-cocycleは中心を既に除いた正しい
linearizationの段階で消えており、projective bundleのBrauer classを零と
根拠なく仮定していない。元稿のtorus説明も、実ruling étale algebra Eに対する
Res_(E/B)GmのtorsorとHilbert90として整合する。この直接vector descentを
本レビューの主証明とする。

## projective no-nameと正確なlevel

W4=Std4のS5作用も一般自由である。
`Y'=(P(V5)×W4)/S5` をP(V5)/S5上で見ると、線形W4からのassociated
rank4 vector bundleなので C(Y')=B(t1,...,t4)。
一方、W4/S5の体はC(e2,e3,e4,e5)で有理4次元。
その上ではY'は、**線形V5表現から降下したrank5 vector bundle**の
projectivizationである。真の線形liftを上で確認したので、generic fiberは
Brauer-twisted形式ではなく普通のP4。ゆえにC(Y')は有理8変数体となる。
従ってBのstable level≤4、F=B(s1,s2)のstable level≤2である。

この計数はprojective P4を用い、そのfiber次元は4。affine V5の余分な一変数を
数えない。FやP4/S5自体の純有理性、最小stable level、多項式不変式環の全生成、
NEW25の方程式は未証明のままである。
