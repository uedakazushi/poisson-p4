# 本質三変数の既約退化五次：全cotangent源と射影次元23上界

2026-11-07、PRIMARY G。新しい非線形計算を行わず、全五変数のcanonical liftから証明する。
命題は一般点の標本・固定H切片・選択したSaito源を仮定しない。
最後の応用で、rank二plane・位数四・squarefree非特別crossratioの
旧concurrent上界31枝を空にする。特殊crossratioのrank二pencil等は残す。

## 定理と適用範囲

射影空間 P4 上の一般階数四 Poisson 構造の実退化五次式を D とする。
三次元の線形部分空間 W⊂C5* があり、D は W の変数だけの既約被約五次式 F であるとする。
S=Sym W とおく。

1. F は exponents (2,2) の自由平面五次であり、全Tjurina数は12である。
2. 任意の線形補空間を固定し、残る二変数を u,v とすると、canonical divergence-zero lift は

       Π=U∧∂u+V∧∂v+B∂u∧∂v,
       B=−(div V)u+(div U)v+R,

   と書ける。U,V は S 上の二次 gradient syzygy の自由基底で、
   R は W の任意の二次式である。vertex変数への隠れた依存はない。
3. この因子型を持つ全点の、PGL商を取る前の射影像の次元は23以下である。
   従って未知G優越成分の一般元にはならない。

自由性と全cotangent表示の証明は、Dが被約で本質三変数なら可約の場合にも成立する。
局所型の有限性と次元23上界の証明には、Fの既約性を用いる。
本定理から、これら全点の既知成分へのE帰属は主張しない。

## 1. canonical lift と Casimir 恒等式

座標を (a0,a1,a2,u,v)、Euler場を E、定数体積形式を Ω とする。
射影二ベクトルを表す二次斉次五変数liftは、E∧X（Xは一次場）を加える自由度を持つ。
divergenceは一次・trace零の場であり、trace零Xについて
div(E∧X)は符号規約を除き5Xなので、divergence零の代表Πが一意に取れる。
scalar Euler場に対応するXの自由度はE∧E=0で消える。
ここではこのcanonical liftだけを使い、任意のEuler代表のrowをlogarithmicとは置かない。

このΠはstrict Poissonである。実際、射影Jacobiは E∧[Π,Π]=0を与え、
divΠ=0とBV恒等式からdiv[Π,Π]=0である。
β=ι_[Π,Π]Ωは重み5の二形式で、ι_Eβ=dβ=0となる。
Cartan公式5β=dι_Eβ+ι_Edβからβ=0であり、[Π,Π]=0が従う。

α=ι_(Π∧Π/2)Ω と定める。これは係数次数四、全Euler重み五の一形式である。
strict JacobiとdivΠ=0のBV恒等式からdα=0である。
実退化式の規約を D=ι_Eα とすれば、Cartan公式は

    5α=dD

となる。奇数次反対称行列のPfaffian恒等式 Π#α=0 により

    Π#dF=0.                                              (1)

を得る。これは既存の `DOUBLE_QUADRIC_SUPPORT.md` 第1節等の標準Casimir恒等式と同じもので、
ここでは必要なdivergence/Jacobi/degreeの根拠を明記した。
canonical liftの発散から射影stabilizerを推論する議論は用いない。

## 2. 全base–base括弧は零になる

Πのbase三行三列をAとする。Aは全五変数の二次斉次係数を持つ反対称行列で、
(1)から A∇F=0である。
Aの三つの係数を符号付きvector b とすれば、b×∇F=0となる。
fraction field上で b=c∇F と書ける。

Fが被約、標数零なので、三つの偏微分F_iの最大公約因子は1である。
実際、共通prime因子はEuler公式5F=Σa_i F_iからFを割る。
Fの各被約因子のgeneric smooth pointではその因子の微分が非零なので、
その因子が全偏微分を割ることはない。
最大公約因子1はvertex二変数を加えた多項式環でも保たれる。
UFDの各primeのvaluationを比較すると、cは分母を持たない多項式でなければならない。
しかし b の次数は二、∇Fの次数は四なので、b=0である。
従ってA=0であり、最初からvertex依存を捨てずに

    Π=U∧∂u+V∧∂v+B∂u∧∂v                            (2)

となる。ここでU,Vはまだ全五変数に依存する二次base場、Bも全二次式である。
(1)はU(F)=V(F)=0を与える。

## 3. vertex零評価から完全なSaito基底を得る

(2)を二乗してEulerを掛けると、fiber成分のEulerは消え、
選んだ順序の符号規約で

    D=−det(E_base,U,V).

これは全五変数上の恒等式である。u=v=0へ評価し、U0,V0と書けば

    det(E_base,U0,V0)=−F.                              (3)

U0,V0はbaseだけの二次場で、Fを消す。
Fが被約なのでSaito determinant criterionは

    Der_S(−log F)=S E_base ⊕ S U0 ⊕ S V0

を与える。EulerのFへの作用は5F、他の二本はFを消すので、
gradient syzygy全加群は正確に

    AR(F)=S U0 ⊕ S V0.                                (4)

である。全基底が二次なので、次数零・一次のgradient syzygyは無い。

全vertex monomialでU,Vを展開すると、その各係数は次数二以下のgradient syzygyである。
正のvertex次数の係数はbase次数一以下なので、(4)により全て零である。
従ってU=U0,V=V0である。
同じ結論は(4)をS[u,v]へ平坦拡大し、全次数二から係数が定数と比較しても得られる。
全fiber依存を初めから除外した選択源の議論ではない。

## 4. Jacobiが残る全係数を決める

U,Vは二次gradient syzygyの基底だから

    [U,V]=αU+βV

でα,βはbaseの一次式である。
(3)にU,VのLie微分を作用させると

    α=div V,   β=−div U.

が従う。例えばU(F)=0とLie微分の体積への作用から、V係数は−div Uとなる。
strict Jacobiの三つ組(a_i,u,v)は正確に

    [U,V]+B_u U+B_v V=0.

U,Vの独立性から B_u=−α、B_v=−βであり、

    B=−(div V)u+(div U)v+R(a0,a1,a2),

Rは任意の二次式となる。
残るJacobiの三つ組は、base括弧零とU,Vのvertex非依存から全て恒等零である。
したがってこれは必要条件だけでなく、全六磁気係数を保持する完全な表示である。
符号は(2)のU∧∂uという規約による。

## 5. Jacobian idealの飽和とTjurina数

Saito基底の次数は1,2,2で、projective exponentsは(2,2)である。
J_F=(F_0,F_1,F_2)の全syzygyは(4)なので、graded resolutionは

    0 → S(−6)² → S(−4)³ → J_F → 0.                  (5)

J_Fはheight二のperfect ideal、S/J_Fは一次元Cohen–Macaulayである。
従ってirrelevant idealに支えられたtorsionは無く、J_Fはsaturatedである。
FがEuler関係でJ_Fに入るため、対応する有限射影scheme Σは全特異点のTjurina schemeである。
Milnor algebraのHilbert numeratorは

    1−3t⁴+2t⁶

で、二階微分の半分は12である。従って length(Σ)=τ(F)=12となる。
このτはMilnor数と区別し、まだ局所型の同次性を仮定していない。

## 6. 既約五次全族の次元は八以下

Fの正規化種数をg、各特異点の枝数をr_pとする。
算術種数は六なので、Milnor–delta公式は

    μ(F)=12−2g−Σ_p(r_p−1)≤12.

各局所でτ_p≤μ_pであり、全τ=12だから、g=0、全r_p=1、全μ=τ=12である。
従って全特異点はunibranchであり、各々でμ_p=τ_pである。
ここで用いる標準定理は、孤立複素平面超曲線germについてのSaitoの
「μ=τなら解析的座標変更でweighted homogeneous」という判定である。
これはfree divisorのdeterminant criterionとは別の依存結果である。

既約weighted homogeneous平面germは、解析的に x^a+y^b、
gcd(a,b)=1、a,b≥2の形になる。
そのMilnor数は(a−1)(b−1)である。
全μ=12なので可能な(a,b)は有限で、複数特異点への分配も有限である。
各型は有限決定的であり、有限jetの座標変換群のorbit条件と特異点の位置のincidenceを使って、
全族は有限個のconstructible equianalytic strataで扱える。
ここで全rational cuspidal quinticの大域分類や、そのPGL軌道の有限性は仮定しない。

(5)をsheaf化し五だけtwistすると

    0 → O_P²(−1)² → O_P²(1)³ → I_Σ(5) → 0.

従ってH¹(I_Σ(5))=0、h⁰(I_Σ(5))=9である。
特異点の位置移動を含むequianalytic接空間は

    H⁰(I_Σ(5))/C F

に入るので、その次元は八以下である。
局所解析型の有限stratificationと合わせ、固定W内の全projective F族の次元は八以下となる。
全ての平面五次形式20次元から、Tjurina十二条件が独立であることを(5)で確認した形である。

## 7. 全Poisson射影像は23次元以下

固定したprojective Fについて、(4)の全frameはGL2の四係数である。
JacobiはBのfiber一次係数を一意に定め、Rは六係数である。
全sourceはGL2×A6で覆われ、tensorの射影scaleを除くと九次元以下となる。
固定projective Fを用いているため、det(frame)が与える実Dのscalarは追加parameterとして数えない。

三次元本質変数Wの選択はGr(3,5)の六次元である。
各Grassmann chartで一つの補空間とframeを固定すれば、定理の全cotangent表示は
その座標での任意のΠを既に含む。
vertexのbase一次shearを別の六parameterとして足してはいけない。
同じshearはU,VのframeとRの変化として、この完全表示に含まれる。
したがって、全sourceの射影像は

    6 + 8 + (4+6−1) = 23

次元以下である。これはPGL商後の次元ではない。
全像の閉包の次元も同じ上界なので、未知G一般成分の25〜27次元という必要条件に反する。

## 8. rank二plane・m4非特別squarefree枝への全適用

既済の `PLANE_REDUCIBLE_CUBIC_M4_SCALAR_BOUND.md` 第1節と独立
`PLANE_REDUCIBLE_CUBIC_M4_SCALAR_REVIEW.md` 第2節を入力とする。
plane上のrankは二、plane沿い実Dの位数は四、leading binary quarticはsquarefree、
classical invariants I,Jは共に非零とする。
全nonlog接続も含むparallel leading tensorの議論は、全concurrent型にも

    G4=L*g4(z,w),

を強制する。ここでLはplane上の非零一次式、g4はconstant squarefree binary quarticである。
五次斉次式のplane ideal I=(z,w)に関する次数から、実D全体は必ず

    D=L*g4(z,w)+g5(z,w)                                (6)

となる。g5は任意のbinary quinticで、未知のnormal係数を固定しない。
Lを三番目の変数とすれば、(6)は本質三変数以下である。
この構造はP1/P2の非線形compatibilityや旧rank3 linear mapのkernelを解くことに依存しない。

G一般元に必要なDの既約性を仮定すると、gcd(g4,g5)=1である。
本質二変数以下の斉次五次はC上で一次因子へ分解するため、Dは本質三変数ちょうどとなる。
従って前節の全G上界23を直接適用できる。

さらにこのsquarefree m4枝は、rank四でD既約という仮定の下では実際に空である。
平面curve (6)のp=[L:z:w]=[1:0:0]では、最低次部g4は四相異なる接線を持つ。
従ってpはordinary quadruple pointで、δ_p=6、μ_p=9である。
既約五次の算術種数六をこの一点が使い切るので、他の特異点は無い。
全τ≤全μ=9である一方、定理はτ=12を強制する。矛盾する。

旧concurrent源の粗上界31を、さらに小行列や二次係数を増やして改善したのではない。
全Dの構造と普通四重点のgenus budgetで、当該全枝を一度に排除した。
同じ適用は旧非特別三型も含むが、その受理済み上界を欠陥なく再計算する必要はない。

## 残余と証明状態

新しい計算証明書やproducerの変更は無い。これは全源に対する幾何・graded moduleの証明である。
独立数学reviewの結果は別ノートに保存する。
二つのSaito定理、有限決定性とequianalytic接空間、平面curveのgenus/δ/Milnor公式は
標準定理として使用し、その適用仮定を上で明記した。

rank二planeのspecial crossratioでleading pencilがrank二になる全枝は、
Dが本質四変数となり得るため(6)へ還元していない。
leading quarticの重根、m2/m3、非被約cubic R、他のrank二/特異/vertical支持、
およびG排除と異なるE全点帰属も未解決のまま保持する。
全成分数や新成分の登録は行わない。
