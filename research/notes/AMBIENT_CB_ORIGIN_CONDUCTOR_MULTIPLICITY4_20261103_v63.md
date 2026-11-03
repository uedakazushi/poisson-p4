# CB原点conductorの重複度4・全initial・depth21

2026-11-03。独立監査PASSの**新しいCB実selfpair Q_b**を使う。
CAのP(a)をCBの二旗証拠として流用しない。
all-nonzero Q_bの追加分類もこのgeneric multiplicity証明には必要ない。
本companion自体の独立監査は別に記録する。

## 1. 独立なCB算術入力

Q_b=P0−b(2C[1]+C[9])は全physical150係数で二つのCB sourceに入る。
独立監査された入力は、両source微分rank29、同じ点へのtangent union rank37、meetingのembedding dimension21、GLとparameterによるimage cone次元20である。
P0の全proper fiber長2の開集合から、generic b≠0の全fiberは二つの被約点となる。
従ってこのgeneric点の全CB completed imageとnormalizationは

    Ahat=T/(I1 intersect I2),       Bi=T/Ii=C[[29]]

であり、第三のsource branchは落としていない。

## 2. intrinsic conductorを使うmultiplicity下界

全射Bi→D_pair=T/(I1+I2)のfiber productがAhatである。
normalizer conductorはker(B1→D_pair) directsum ker(B2→D_pair)なので、
intrinsicなAhat/conductorはD_pairそのものとなる。
従って実meetingのedim21は、別のquadratic R29 chartとの座標同一視を仮定せず、actual conductorのedim21を与える。

受理済みupper theoremによりconductorの次元は20である。
実Q_bのG-saturation Zも20次元でP0を通るため、conductorのsupport componentとなる。
もしZ上でgeneric conductor長が1なら、regular locusの稠密開部分ではedim20となる。
しかしconstructibleで稠密なG×b族はその開集合に会い、そこでも受理済みedim21が成り立つので矛盾する。
従ってZのgeneric長は2以上である。これは「二旗だから非被約」という推論ではない。

generic Q_bの二つのnormalization点は同じR値と異なるt値を持つ。
B=R[t]/(t²−nu)でtが全Bを生成するので、nu=(t1−t2)²/4は非零である。
従ってZはnu=0のcritical W、V componentと異なる。
両flag arcはpropernessとP0の一台fiberによってP0の完成近傍へ延び、そのclosureは原点を通る。

W、Vは既に証明されたgeneric Ac(8)型から各々smooth20・generic長1である。
associativity formulaによりe(conductor)≥1+1+2=4となる。
completionでsupportが分岐しても、この下界は弱くならない。

## 3. 全initialとCM性

先行theoremは全minimal initial idealに

    I0=(WU,WV,U*(aU²+bUV+cV²)),       a!=0

が含まれ、I0の商がCM20・multiplicity4であることを証明した。
upper4と上のlower4が一致するのでe(conductor)=4である。
全射S_gr/I0→gr(conductor)のkernelが非零なら、CM・unmixedなsourceの非零submoduleとして正の20次元multiplicityを持つ。
これはmultiplicityの等号に反する。従って **in(I')=I0** である。
b,cを未計算のまま零と置く操作はない。

Hilbert seriesは(1+T)²/(1−T)^20であり、conductorはCM20、edim22、multiplicity4、type2である。
minimal support primeは正確に三つ、generic長は(1,1,2)、embedded associated primeはない。
generic Q_bの第三component上ではconductor完成環はC[[g20,eps]]/(eps²)となる。

## 4. 全actual odd idealとCB imageの不変量

三つのinitial generatorのactual liftsがI'を全て生成することは、complete local divisionから従う。
七actual regular odd functionを戻すとmu(I)=10であり、Iはactual coefficientの次の十行で生成される。

    18,19,28,29,60,67,72,46,95,63.

これは全analytic normal formや混合cubic係数を計算したという意味ではない。
正則R29上でconductor R/IがCM20なのでdepth I=21、A=R directsum I*tからdepth A=21となる。
mu(I)=10によりedim A=39であり、odd productsがm_R³に入るためmultiplicityはR-rank2である。
coneではdim29/edim39/depth21/e2、projectiveではdim28/edim38/depth20/e2となる。
原点のactual CB completionはregular normalizerのsubringとしてdomainのままである。
conductorの非被約構造をCB imageのnilpotentとは解釈しない。

CAの独立な対応theoremと結果は一致するが、selfpairの入力は異なる。
full binary cubic係数・higher analytic同型型、actual CA/CB intersection D、full Poisson schemeはここでは決定しない。
