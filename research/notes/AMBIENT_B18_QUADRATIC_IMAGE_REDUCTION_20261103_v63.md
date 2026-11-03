# generic B18におけるactual CA/CB交叉の共通rank2 image-kernel還元

2026年11月3日。D=CA∩CBは**被約成分どうしの実際のscheme intersection**である。
本稿はgeneric B18での全nilradicalやprimary decompositionを決めるものではない。
代わりにactual Dを、同じregular28変数底上の二つのrank2有限algebraの**image kernelsの和**へ正確に還元する。
normalizerどうしのtensor productでactual Dを代用しない。

## 結果

B18の非空稠密PGL5軌道の代表点を

    P0=∂y∧D0,
    D0(h)=x², D0(x)=0, D0(z)=hx, D0(w)=h²+z²

とする。この点はB17外にある。

1. 全proper CA/CB/N sourceのP0でのfiberは、いずれも一つの台を持つ長さ2のschemeである。
   各sourceはtargetの近傍で有限normalizationとなる。
2. **共通のactual target座標28個**からなる正則底R=C[[r1,...,r28]]を具体的に選べる。
   CAとCBの完成normalizerは、ともにR上**finite free rank2**である。
3. 同じtarget regular ring T=R[[z1,...,z97]]について、全actual完成交叉は正確に

       Dhat = T/(ker phi_CA + ker phi_CB)

   と書ける。各phi_iは明示的な有限source circuitから定義され、targetをR-rank2のnormalizerへ送る。
   各kernelは、下記の**一変数quadratic algebraへのimage kernel**、すなわち一変数eliminationである。

したがって元の二つの28次元smooth sourceを丸ごとtensorして得るschemeを用いず、
二つのdegree2 image-kernel計算が残る有限なactual問題となった。
quadratic係数やodd idealの全展開、二つのkernelの和のprimary構造は未計算である。

## 1. 全geometric旗と全scheme fiber

既存のproper模型により、CAのrank≤2はg=0またはf=0、CBではassociated ratioのendpoint二つに限られる。
P0の定covector Casimirはxの一つだけである。
f=0のCA endpoint、およびa=0のCB endpointは少なくとも二つの独立Casimirを持つため、この点には現れない。
従ってCAではg=0、CBではb=0側を扱えばよい。

この二つのendpointはcanonical代表で定vector factor Lを持つ。
P0についてE∧v∧P0=0を満たす定vectorの直線はL=C∂yだけである。
x=0に制限するとD0(w)=h²+z²が非零なので、商vはw方向に限られ、D0(h)=x²でその余りも0となる。
この条件はLに線形でrank4だから、fiber上のLはschemeとしてもreducedな一点である。

残るmarked covector ellはD0(ell)=m*ellを満たすlinear Darboux formである。
ellにw成分があれば、w²・w係数の比較でm=0となり、次にD0(ell)=0から矛盾する。
ell=A h+B x+C zならD0(ell)=x(Ax+Ch)なので、可除性からellはxに比例するしかない。
従ってmarked lineはxで一意である。
そのcofactor m=0からCA/CBの両方でf_k=0となり、D0(k)はx²の定数倍である必要がある。
これよりW=span(x,h)も一意に定まる。以上はchart内だけのuniquenessではない。

この唯一の旗の近傍で、残る全五旗変数を

    ell=x+a*h+b*z+c*w,
    k=h+d*z+e*w

とする。L固定後のF(1,2;4)の全chartである。
有限coordinate matrixを作用し、canonical150係数が実際のfixed-flag source空間へ入る全linear constraintsを作った。
QQ Groebner basisは正確に

    CA: (b,c,d+2a,e,a²),
    CB: (b,c,d−a/2,e,a²),
    N : (b,c,d+2a,e,a²).

rank2条件は各sourceのf-unit chartでgまたはratioをschemeとして0にするので、
このfiber計算でCBの非線形ratio条件を落としてはいない。
各fixed-flag tensor embeddingはinjectiveで、fiber係数はtargetから一意に戻る。
従ってこれが全proper source fiberで、長さは2である。
propernessで有限化し、既存のbirational smooth sourceからfinite normalizationを得る。
一つのnormalizer点しかないので完成targetはその正則完成環へ注入する。

## 2. 二つの全normalizer chartを有限に作る

sourceの座標順序(h_N,k_N,y,z,w)=(x,h,y,z,w)を使う。
CAの20個のfixed-flag基底はfの三binary方向、g h_N²の一方向、Hの16方向である。
CBではf三方向、g h_N²の接方向、H16を使うが、全finite族では必ずg=mu*fを課す。

中心のHは−h_N*k_N*w+k_N²*z+z³/3である。
CA、CBのvertical旗方向を与えるLie行列はそれぞれ

    MA=E_xh−2E_hz,
    MB=E_xh+(1/2)E_hz.

いずれも三乗が0で、exp(pM)=I+pM+p²M²/2は有限polynomialである。
source fiberを−p*Lie_M(P0)だけ補正し、保存した九fiber基底を独立に動かすと、
全十変数(p,k0,...,k8)のfinite core mapを得る。
CBではg=−alpha3*fを用い、mu*fの全cross termsを含めた。
全150出力を証明書へ保存した。

さらにP0のGL orbitから十九の独立方向を選ぶ。
それらのelementary GL行列の有限積H19をcore mapに作用させたものが全source chartである。
source imageの一次rankは28で、pのtarget微分は0である。
一方pの旗微分は非零であり、旗Pluecker係数を加えた全29 source微分はrank29である。
よってH19×core10はsmooth cone-normalizerへのetale chartで、全source completionを与える。
この段階でsource image rankをtarget tangent dimensionと取り違えていない。

## 3. 共通target底と二次有限性

raw150順序で選んだ共通28行は

    5,15,16,20,22,23,27,31,35,50,65,105,106,108,
    109,110,112,113,117,118,119,120,121,123,125,127,132,140.

源のp以外の28方向に関するminorはCAで−2304/25、CBで−1536/25である。
従って両方で同じ28実target座標を固定して、残る28 source変数をimplicitに消去できる。

core mapのp²係数Wについて、このimplicit eliminationによる正確な補正は

    W' = W−J*(J_selected)^−1*W_selected.

29個目のtarget座標にraw row61を選ぶと、W'のその係数は

    CA: −4,   CB: 3/8

で、両方とも非零である。
従って共通29 cone target座標を0としたnormalizer fiberは正確にC[[p]]/(p²)となる。
高次termは単元倍のp²へWeierstrass化され、長さを増やさない。

Rcone=C[[29 actual target coordinates]]とすると、complete Nakayamaで各normalizerはRcone上有限である。
normalizerは正則29次元で、これら29個はparameter idealをなす。
CM性とregular底によりfinite moduleはfree、closed fiber長が2なのでrankは2である。
これは単に全source fiberが長さ2ということからfinite flatを推論する議論ではなく、
**共通のactual parameter idealに関する長さ2**を別に確認している。

row106は選択28行に含まれ、中心で1である。u106=1の射影chartを取ると、
同じ共通底の一変数を除けてR=C[[r28]]となり、両normalizerは依然free rank2である。
canonical150空間の24本のlinear divergence constraintsを除き、選択29行を126個のraw座標基底へ延長した。
その全constant reconstructionを保存したため、共通target ringは射影で正確にT=R[[z97]]となる。

## 4. actual image idealsへの正確な一変数還元

各normalizer BiはR上free rank2で、閉fiber basisが1,p_iであるため、同じ二つがR-basisとなる。
従って一意な係数tau_i,nu_i∈Rがあり

    Bi=R[[p_i]]/(p_i²−tau_i*p_i−nu_i).

各残りの97 actual target座標のBiでの像は一意に

    phi_i(z_j)=a_ij(R)+b_ij(R)*p_i

となる。これらの全係数は次の**有限な代数的presentation**で指定される。
上のH19/core10有限source circuitへ29 common target方程式を課し、
28個のunit-linear pivotをHensel消去して残る一つのWeierstrass二次式を取る。
全target出力をそのrank2 basisで一意にdivideする。
係数はalgebraic formal power seriesとして指定され、未検査のjet切断ではない。
展開表を計算済みとは言わない。

したがって

    ker phi_i = (p_i²−tau_i*p_i−nu_i,
                 z_j−a_ij−b_ij*p_i for all j) intersect T.

これは**一変数の二次elimination**である。
各contractionを別々に求め、その**和**を取るのがactual Dである。
二つのnormalizer graphを先に足してから同時にeliminateしたり、Biどうしのtensorをactual Dと同一視しない。

同値にIi=(b_i1,...,b_i97)⊂Rとすればactual component imageはR⊕Ii*p_iである。
その全kernelはIiのfirst syzygiesと、y_j=z_j−a_ijに対するpair products

    y_j*y_k−tau_i*b_ij*y_k−nu_i*b_ij*b_ik

で尽くされる。これらで全polynomialをR+ΣR*y_jへ戻し、残るkernelはrank2 basisの独立性からsyzygiesだけになる。
従って残る計算は二つのfinite ideal Iiとその実target埋込みの比較であり、
normalizerの余分なnilpotentsをactual intersectionへ混ぜていない。

## 証明書と未解決部分

* `code/ambient_B18_quadratic_image_reduction_20261103_v63.py`
* `reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json`

全QQ計算3.21秒。全flag constraints/GB、両finite150 core maps、全29 augmented source微分、
共通28 minor、全150 corrected quadratic係数、126座標へのconstant reconstructionを保存した。
P0のprojective orbit次元18と既存B18次元18から、これはgeneric B18に適用される。

actual Dのgeneric B18での全completed ring、nilradicalの冪零次数・annihilator・embedded primesは未確定である。
二つのquadratic-image係数idealを実際に求めることが次のbounded課題であり、今回のreductionを解決済みprimary typeと呼ばない。
