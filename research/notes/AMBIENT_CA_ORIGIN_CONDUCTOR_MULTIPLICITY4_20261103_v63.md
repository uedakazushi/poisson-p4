# CA原点conductor：重複度4、全initial、CM性とactual imageのdepth21

2026-11-03。**CAだけ**の定理である。
両CA/CBの受理済みupper20/e≤4と、独立監査済みCAの実selfpair族を結合する。
二つの異なる旗だけからgeneric conductorの非被約性を推論しない。
CB、CA/CBのactual intersection D、full Poisson schemeへ移す主張はない。

## 1. actual selfpairのmeetingがintrinsic conductorとなる

実selfpair族P(a)はP0を通り、GLとparameterによるsaturation Zはcone20次元である。
P0を含むproper有限fiber長≤2の開集合上で、generic a≠0の全fiberは二つの被約点である。
両sourceはsmooth29かつactual source微分rank29である。
受理済みの結論として

    Ahat=T/(I1 intersect I2),
    Bhat=B1 directsum B2,       Bi=T/Ii=C[[29]]

が全completed CA imageと全normalizationを与える。
whole Aのembedding dimensionを接空間の和37と置く仮定は使わない。

D_pair=T/(I1+I2)と置けば、全射Bi→D_pairとexact sequence

    0 -> Ahat -> B1 directsum B2 -> D_pair -> 0

がある。Ahatはこの二つの全射のfiber productである。
normalizer conductorはker(B1→D_pair) directsum ker(B2→D_pair)であり、
従ってintrinsic conductor quotientは **Ahat/conductor = D_pair** となる。
これはcanonicalな同定であり、以前のquadratic base R29と別chartの座標を同一視する追加仮定を必要としない。

## 2. generic長が2以上である理由

D_pairのembedding dimensionは、全QQ(a) tangent計算によって21と証明されている。
一方、先行upper theoremによりconductor全体の局所次元は20である。
Zはその20次元のirreducible saturationなので、conductorのsupport componentとなる。

もしconductorがZのgeneric pointで被約なら、他成分とnilpotentのsupportを除いたZariski openではregularとなり、embedding dimensionは20となる。
finite type・characteristic zeroのregular openであり、constructibleで稠密な実G×a族はそこに会う。
しかしこの族のgeneric aと全G移動でembedding dimensionは21であるから矛盾する。
従ってZ上のgeneric Artin長は少なくとも2である。
二つの旗が存在するというだけの推論ではなく、**実edim21・次元20・族の稠密性**を使う。

critical W、Vの二つのsmooth20 support branchはnu=0に含まれ、generic長はそれぞれ1である。
selfpairの二つのsource点は、P0の完成近傍では同じR値を持つ異なるt値を与える。
B=R[t]/(t²−nu)だから、二点が異なることとtがR上の全Bを生成することから

    nu=(t1−t2)²/4 != 0

がgeneric selfpair arc上で成り立つ。
従ってZはcritical W、Vとは異なる第三のcomponentである。
二つのflag arcはa=0でP0へ合流するので、そのclosureは原点を通る。

## 3. e=4が全initialを固定する

associativity formulaを原点のcompleted conductorへ適用する。
WとVのsmooth branchの寄与は1ずつで、Zからは少なくとも2が寄与する。
Zがcompletionで分岐してもこの下界は弱くならない。
従ってe(C)≥1+1+2=4である。
先行upper theoremのe(C)≤4と合わせ、**e(C)=4** となる。

先行theoremは、minimal embedding S22でgenuine initial inclusion

    I0=(WU,WV,U*(aU²+bUV+cV²)) subset in(I'),      a!=0

を与える。I0の商はCM20次元・重複度4である。
全射S_gr/I0→gr(C)のkernelが非零なら、CMかつunmixedなsourceの非零submoduleとして20次元を持ち、正のmultiplicityを持つ。
これは両商のmultiplicity4の等号と矛盾する。
従ってkernelは零で、**in(I')=I0** である。
この議論でb,cや未知higher係数を零に置かない。

Hilbert seriesは

    H_gr(C)=(1+T)²/(1−T)^20

である。従ってCはCM、dim20、edim22、multiplicity4、CM type2である。
minimal support primeは正確に三つ、generic長は(1,1,2)であり、Cにembedded associated primeはない。
Zのreduced component自身の原点multiplicityも1に確定する。

## 4. actual IとCA imageへの帰結

I0は二つの二次式と一つの三次式で生成される。
initial equalityとcomplete local divisionにより、そのactual lift三つがI'を全て生成する。
元の七regular actual odd functionを戻すと **mu(I)=10** であり、次の十actual odd coefficientがIを生成する。

    rows 18,19,28,29,60,67,72,46,95,63.

これは全analytic normal formやbinary cubic係数を求めたという主張ではない。
選択十functionは既存finite trace recipeで定義された実際のfunctionであり、抽象的な未指定generatorではない。

正則R29上でC=R/IはCM20なので、depth I=21である。
actual A=R directsum I*tから **depth A=21** となる。
mu(I)=10、odd productsがm_R³に入ることから、Aのembedding dimensionは29+10=39、multiplicityはR-rank2である。
従ってconeではdim29/edim39/depth21/e2、projectiveではdim28/edim38/depth20/e2となる。
原点のAはregular normalizer Bのsubringであり、domainであることは変わらない。
Cの非被約構造をAのnilpotentと取り違えない。

一般selfpair点ではW/V componentを避け、Cはsmooth20 support上のgeneric長2のCM schemeとなる。
この開集合の完成環はC[[g20,eps]]/(eps²)であり、二smooth29 CA branchはこのordinary double20 meetingを持つ。
これは個別CAのselfintersectionの話で、CAとCBのintersection Dの型ではない。

CBについては現在、同じupper20/e≤4までを別theoremとして保持する。
新しいCB selfpairが構成されても、全physical source・proper fiber・edim21の独立監査前にはこの等号の議論をCBへ移さない。
