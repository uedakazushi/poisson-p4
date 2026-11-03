# CA・CB原点conductorのbinary cubic initial：独立認識

2026-11-03。元のactual physical selfpair arcをdecodeして計算し、既知fixed-basisの短縮式を入力にしない。
既に受理済みの全initial・CM性・三support branch・generic長(1,1,2)を使う。
ここで確定するのはassociated graded conductorの三次式であり、全analytic normal formではない。
凍結済みの「混合係数未決」と記したhistorical sourceやreportは変更しない。

## 全29式での独立tangent計算

CAでは凍結されたP(a)の全150係数を、degree/rational-coefficient形式から直接decodeする。
CBでは凍結されたQ_bの全150多項式を直接読む。
各parameterを0で微分し、共通actual R29座標に制限する。

原点upper theoremのK29×22は、十九G方向とW,U,V方向を順に並べたminimal tangent kernelである。
Kの一つのsquare minorだけを使うのではなく、全29式K*c=rを有理数Gaussian eliminationで直接解き、全29残差が厳密に零であることを検査する。
得られたWUV成分は

    CA: (0,−5/2,−15/2),
    CB: (0,−9/4,−1/4)

である。従ってself supportのtransverse tangent lineはそれぞれ

    CA: U−V/3=0,       CB: U−9V=0.

G19とこのarc微分を合わせたrankは、原点そのものでも20である。
そのQQ minorも保存するため、self supportの原点での滑らかさをgeneric family rankから外挿していない。
critical W、Vのsmooth supportと合わせ、三つのtransverse tangent lineは互いに異なる。

## 長の初期形への配分

intrinsic conductorはG不変である。
全150係数のsmooth ambientで、stabilizerに横断する十九formal group方向と補sliceを取れば、action mapはIFTによりformal isomorphismとなる。
G不変idealのpullbackはgroup変数に依存しないので、conductorを十九smooth因子とcurve coreに分けられる。
R29自体のG不変性や任意image/base-change交換は仮定しない。

curve coreはCMであり、三つのreduced support branchは原点で滑らか、接線は互いに異なる。
generic primary長は順に1,1,2である。
prime filtrationとRees specializationにより、tangent coneのfundamental cycleは

    [W-axis]+[V-axis]+2[selfline]

となる。各smooth curve branchはDVRであり、誘導されるgood filtrationはgeneric tangent line上でrank1を持つ。
したがって各branchのgeneric長が、その異なるtangent lineのcycle係数として現れる。
nilpotentを含むというだけで長2を別の方向へ割り当てたり、総multiplicity4だけで重複因子を推測したりしない。

受理済みの全initialはCM・unmixedな

    (WU,WV,U*(aU²+bUV+cV²)),       a!=0

であり、追加のinitial componentはない。
従ってplane cubicにはV-axisの単純因子Uと、selflineの二重因子が必要である。
保存したactual U³ classが全体の定数を固定するので、正確に

    CA: f3=(12/125)*U*(U−V/3)²,
    CB: f3=(−4/243)*U*(U−9V)²

となる。
展開した(U³,U²V,UV²)係数は

    CA: (12/125,−8/125,4/375),
    CB: (−4/243,8/27,−4/3)

である。

この認識は全initial idealを具体化する。
元のactual conductor idealのhigher項、normalization coverの全analytic同型型、CAとCBのrelative embeddingやactual intersection Dは決定しない。
