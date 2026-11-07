# 本質三変数の既約五次：全磁気項除去と有限23軌道への入口

2026-11-07、PRIMARY Eの全源認識。
`ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md` の全canonical cotangent定理を入力とする。
本稿は任意の固定curveで全六磁気係数を実projective shearによって消し、
本質三変数の既約被約五次を持つ全rank四Poisson支持を有限個の23次元軌道閉包へ還元する。
曲線の具体的列挙、各軌道閉包の成分性又は既知成分への帰属は未解決である。
これは全P4の軌道分類を新たな完了条件とするものではなく、この残E支持の有限全源入口である。
新CAS・jet・点列・minor計算は行わない。

## 1. 固定Fの全源と、実shearの作用

S=C[a0,a1,a2]、Fは既約被約五次、Eはbase Euler場とする。
Fはfree exponents (2,2)であり、そのgradient syzygyの二次自由基底をU,Vとする。
前稿の定理により、任意のvertex座標u,vを固定して、全Poisson源は

    Π_R=U∧∂u+V∧∂v+
        (−αu−βv+R)∂u∧∂v,
    α=div V,  β=−div U,  R∈S2

の形で、U,Vの全GL2 frameと任意の六係数Rに尽くされる。
base–base括弧は零であり、U,Vはvertex変数に依存しない。

l,m∈S1を任意に取り、実際の線形座標変更

    u'=u+l(a),   v'=v+m(a),   a'=a

を行う。この5×5行列はunipotent、determinant1で、全parameterで可逆である。
全括弧を変換するとanchor U,Vは変わらず、磁気項は正確に

    R'=R+αl+βm−U(m)+V(l)
      =R+div(lV−mU).                                 (1)

となる。base–base括弧が零なので{l,m}=0であり、これは一次近似ではなく全tensorの恒等式である。
従って以下の線形写像が可逆なら、全Rを一度の実shearで消せる。

    δ:S1⊕S1 → S2,   δ(l,m)=div(lV−mU).              (2)

双方は六次元である。

## 2. kernelからprojective closed一形式を作る

δ(l,m)=0とし、W=lV−mUと置く。
Wは係数次数三のbase場で、W(F)=0、div W=0である。
baseの定数体積形式をΩとし、θ=ι_WΩを二形式とする。
θのEuler重みは五であり、Fの次数も五なのでθ/Fの重みは零である。

div W=0からdθ=0、W(F)=0からdF∧θ=0となる。
よってd(θ/F)=0であり、Cartan公式から

    ω=ι_E(θ/F)

はclosedである。ωはEuler-horizontalかつ重み零なので、P2上のrational一形式に降下する。
その分母はF一回だけであり、polesはFに沿うsimple polesだけである。

Fのsmooth locusで局所座標を(z,t)、F=unit*zと取る。
simple-pole一形式をω=(A dz+B dt)/zと書けば、closednessの二重極部分からBはzで割れ、
残る単純極部分からA|z=0はtに依存しない。
従ってresidueは各smooth branch上の定数である。
Fは既約なので、そのsmooth locus全体で同じ定数cとなる。

一般の直線ℓはFの特異点を避け、五つのsmooth点で横断する。
ω|ℓはP1上のrational一形式で、他のpolesを持たず、五つのresidueは全てcである。
P1のresidue定理から5c=0、従ってc=0となる。
よってωはFのgeneric pointでpoleを持たない。
他のdivisorial polesも無いため、smooth P2の局所自由な一形式層のcodimension二延長により、
ωは大域正則一形式である。H0(P2,Ω1)=0からω=0となる。

ι_Eι_WΩ=0は、EとWがfraction field上で平行であることを意味する。
W(F)=0、E(F)=5F≠0なのでW=0である。
U,Vは自由基底だからl=m=0となる。
従ってδはinjective、六次元同士なのでisomorphismである。

この議論ではlogarithmic comparison theorem、analytic normal form、Fの特異点の分類は使わない。
既約性と被約性、canonical全源、P1上のresidue定理だけを使う。

## 3. 全六磁気項を一度の実projective shearで消す

任意のR∈S2に対し、一意の(l,m)でδ(l,m)=−Rとできる。
(1)の実unipotent shearにより、Π_RはΠ_0へ移る。
磁気零・非零を別のgeneric点とみなしたり、磁気零を単なる軌道極限として扱う必要はない。
全六係数で同じ可逆写像を使う完全なinverseである。

Fを動かす場合、AR(F)_2は全free (2,2) locus上でrank二のkernel vector bundleとなる。
各local frameにおけるδのdeterminantは、上の証明により全閉点で非零である。
従ってlocal inverseはそのframe chartの全域で正則に取れ、未知のspecial Fを除くgeneric inverseではない。
frame chartの補集合は別のframe chartで覆われる。
vertexのcomplement選択も前稿の全canonical表示で既に保持している。

固定Fのframe GL2はvertexの実GL2変更によって互いに移る。
これと全磁気除去を合わせ、固定projective Fに属する全Poisson tensorは一つのPGL5軌道に入る。

## 4. admissible plane quinticのPGL3軌道は有限個

前稿は、Fが有理、全特異点unibranch、全μ=τ=12であることを示した。
Saitoのμ=τ判定から全局所型はweighted homogeneousであり、
irreducible weighted homogeneous plane germ x^a+y^bで(a−1)(b−1)≤12を満たす型は有限である。
有限決定性から、全Fは有限個のconstructible equianalytic strataに分かれる。
前稿の飽和Jacobian resolutionは、特異点の位置移動込みで各stratumの次元を八以下にする。

一方、Saito自由基底の次数は1,2,2なので

    Der(−log F)_1=C E.

従ってPGL3内のFのprojective stabilizerのLie algebraは零である。
標数零でstabilizerは零次元のalgebraic group、ゆえ有限であり、各PGL3軌道の次元は八である。

各stratumはPGL3安定である。
その中の八次元orbitは、八次元以下のstratumのあるirreducible componentでdense openとなる。
同じcomponentに二つの異なるdense open orbitは存在せず、
その境界は八次元未満なので他の八次元orbitを含まない。
finite-type constructible strataのirreducible componentsは有限である。
よって全admissible FのPGL3軌道は有限個となる。

この有限性は曲線の具体的なnormal formリストを既に得たという意味ではない。
また、τ=12だけで任意の曲線がfree (2,2)になるとは用いない。
定理の方向は、Poissonからfree (2,2)を得て、そのlocusにのみ上の議論を適用するものである。

## 5. 対応する全PGL5軌道の次元は23

Fの本質変数部分空間Wは内在的に一意である。
Fが三変数に本質的なので、projective Dを保つPGL5変換はWを保つ。
そのbaseへの誘導PGL3作用のstabilizerは有限である。
従って連結成分内では、全行列のscalarを除いて

    a'=a,    (u',v')^T=M(u,v)^T+C a,

M∈GL2、C∈Mat(2,3)という形に限る。
これはGL2の四方向とshear六方向からなる十次元群である。

Π_0のanchor U,Vは独立なので、これをprojectively保つMはscalar λI2に限る。
shearの磁気変化は可逆なδだから、Π_0を保つにはC=0である。
逆にa'=a、(u',v')=λ(u,v)はΠ_0をλ倍するのでprojective stabilizerとなる。
よって全projective stabilizerの次元は正確に一、identity componentはこのC*である。
従ってactual PGL5 orbitの次元は24−1=23である。

全admissible plane Fは有限PGL3軌道で、各固定Fの全frame・全磁気sourceは一軌道である。
本質三変数部分空間を動かすPGL5作用も含めると、対象の全Poisson支持は
有限個の23次元PGL5軌道の和であり、その有限個の閉包の和に入る。
異なるplane curveのPGL3型は、全Dとその本質変数からPGL5でも区別される。

各軌道閉包は23次元の閉既約Poisson族であるが、より大きい既知又は未知の成分の境界かもしれない。
ここから全空間成分性や成分数の増加を推論しない。

## 次の有限E義務

次は、free (2,2)であるirreducible rational cuspidal quinticの具体的PGL3代表を
一次資料と照合して有限リストにし、各代表の全Poisson軌道閉包を既知成分へ認識するか、成分性を証明する。
Namba、Flenner–Zaidenberg、Dimca–Sticlaruの関連分類は未照合の文献候補であり、
本稿の有限性・磁気除去の証明の入力にはしていない。
未照合のリストを既済分類として引用しない。
既済のGa9 linear-wは退化五次が可約なので、この既約Fの全磁気除去を適用しない。
