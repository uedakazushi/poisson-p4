# Scalar entryの全45一次spanとFitting supportの有限証明書

2026-11-03。5×670 moduleの**個々のscalar entries**を結合すると、全45個のformal parameter一次式を得られることをgeneric betaで証明する。これはvector rowを結合する前の障害とは異なり、**parameterがAをannihilateするという意味ではない**。Fitt4がorigin idealであることのconstructive certificateである。

## 1. generic構造上界179

formal base R'=K[n42,S1,S2,H]のparameterをx_aとする。元qへの代入は

    q_i=P_i(h)+L_i(x), degree_h(P_i)≤4,

であり、L_iはn,S1,S2に一次、Hを含まない。元134式はqに二次なので、parameter二次部分は134本のquadratic form Q_i(L)で、h次数0にだけ現れる。monic relation g(h)−Hで割る際、h≥5の項からHが生まれる。

従ってh^j r_i(j=0..4)のtarget h^k成分の**Hを含まない二次部分**は、厳密にdelta_(j,k)*Q_i(L)である。その他の二次部分は全てH*x_aのspanに入る。よって全3350scalar entryのquadratic span次元は、generic betaでも

    ≤134+45=179

である。これは凍結点のrankからの上界外挿ではなく、全sourceのquadratic形式とmonic divisionの構造による。producerは元134からQ_i(L)を別に構築し、全3350entryのこのpatternをQQ係数で検査する。

## 2. 全1080monomial universeと非零minor

45linear+1035quadraticの**全1080列**を使う。凍結点で零になる六つのquadratic monomialsも列を残し、一般betaでも零だとは仮定しない。237811非零係数をmod32003へ落とし、scalar行列3350×1080の224次square minorを選んだ。その元の行列でdeterminantは**21497≠0**となる。全行列入力、selector、monomial labelsを保存する。

従ってgeneric QQ rankは≥224である。一方、quadratic span≤179とlinear空間次元45から全rank≤224なので、両方一致する。linearに落ちるrow spanの次元は224−179=45であり、全parameter一次式の空間そのものである。選択column basisは実際に全45linear列を含む。

generic K係数は保存224minorを反転して各linear targetを解く有限recipeで指定できる。全rank上界224が証明済みなので、選択列で解いた係数は残る全1080列でも一致する。巨大なQQ inverseをserializeしたとの主張ではない。

最初のprobeは観測1074列の全column rankを期待したが、実rank224であった。全列を残したまま、上のsource構造を用いて正しい上界179と証明を得た。より高次数の探索はしていない。

## 3. Fitting idealsと未知support

正確な5×670 module presentationでscalar entriesが生成するidealはFitt4(A)である。全entriesに定数項はなく、全linearsを含むから

    Fitt4(A)=m=(n42,S1,S2,H)

を得る。これは「scalar entriesを互いに結合する」操作であり、同じtarget成分を持つvector module rowの結合ではない。したがって前のvector-row quadratic rank670、purelinear vector span0と矛盾しない。m⊂Ann(A)も従わない。

さらに既知のscheme splitとnormal閉埋込みより、R'上でA≅(R'/m)^5⊕R'/I_Bとなる。supportsは互いに交わらない。各点でFitting idealsを計算すれば

    Fitt_j(A)=m^(5-j), j=1..4,
    Fitt0(A)=m^5 intersect I_B.

従ってrank≤3 locusはoriginだけで、未知Bの上ではmatrix rank4である。4×4 minorを反転するchartではSchur complementが一行のscalar idealとしてBを正確に表す。未知I_Bの空性はこの議論から従わず、maximal5minorsによるI_B=Fitt0:m^infinityがなお具体的な課題である。

## 範囲と保存

証明は自由formal45 baseであり、非線形S(n),H(n)へのbase change後の定数spanを主張しない。ただしabstractなFitting formulasはnormal42 baseでも同じ直接積から成立する。

`code/s61_affine_normal_scalar_fitting_20261103_v63.{py,cpp}` と `reports/S61_AFFINE_NORMAL_SCALAR_FITTING_20261103_v63/`。一つのbounded60秒・2.2GB rank probeと非零minorの抽出だけで、minor集合列挙・高次数Macaulay・GBは行っていない。
