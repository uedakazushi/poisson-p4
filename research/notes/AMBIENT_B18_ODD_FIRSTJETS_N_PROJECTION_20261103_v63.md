# actual odd係数の一次項と、共通R29へのN像のpinch point

2026-11-03。generic B18点でのactual CA/CB/N imageを調べる。
本稿では、独立な抽象odd ideal型から一歩進め、保存された**共通の実target座標**で
odd係数の一次項を求める。さらにNの共通R29への形式射影像を全次数で同定する。
actual CA/CB交叉Dの全完成環はまだ求めていない。

## 1. common actual coordinatesとsource

既存のCA/CB finite core mapをU(p,k0,...,k8)とし、その外側に19個のGL作用を加える。
U_p(0)=0で、28個の他source方向のtarget微分をJとする。
全150係数を保存し、JのGL19列を元の物理Lie式から作り直して照合する。
sourceは全29次元normalizerへのetale chartであることを、先の独立受理済み証明から使う。

共通28実target行をY、29個目をraw row61とする。
J_Yは両族で可逆である。

## 2. odd係数の真の一次項を有限jetから取り出す

sourceの残る28変数をsとし、展開を

    F(p,s)=F0+J*s+W*p²+Z*p³+P_s*p*s+高次

と書く。P_sのGL19列はU_p(0)=0から零であり、残る9列はU_pkから直接得る。
Yをactual base座標としてsを形式陰関数で解く。
Y=0での二次補正はs20=−J_Y^(-1)W_Yである。
射影rem(V)=V−J J_Y^(-1)V_Yとおくと

    Wc=rem(W),
    Zc=rem(Z+P_s*s20),
    E=rem(P_s*J_Y^(-1))

である。w=Wc_61はCAで−4、CBで3/8となる。
さらにq(V)=V−Wc*V_61/wとおく。
eta=row61−L*Y、L=J_61 J_Y^(-1)というadjusted29番目座標を取れば、
全targetのodd係数の一次項は

    b_linear(Y,eta)=q(E)*Y + q(Zc)*eta/w

である。raw共通R29順序(Y,row61)では、Y列をq(E)−q(Zc)*L/wへ戻す。
各式はfull sourceの必要な次数係数から導いた正確な一次項であり、
任意のjet切断をfull odd idealと同一視したものではない。

理由は、rank2基底1,pでeven部分をsubtractするとき、二次target項のeven寄与は
Wc/w倍の29番目座標であることによる。
pのtraceを引いてt=p−Tr(p)/2へ移しても、b(0)=0なのでodd係数の一次項は変わらない。
座標pの二次変更から生じるWc倍もqで消える。

両族ともY28方向のodd係数微分はrank6、29方向全部ではrank7となる。
従って実odd idealには七つの独立な一次部分があり、そのうち一つはbranch方向を含む。
これはfull odd idealが七生成という意味ではなく、minimal生成数10以上とは両立する。

common R29でのbranch discriminantの微分は

    dnu=(drow61−L*dY)/w.

Lの非零成分は

    CA: L15=−2/3, L22=−1/2, L27=2/3,
    CB: L15=−4,   L22=−1/2, L27=4.

二つのdnuは線形独立である。nuがR/I上の非零因子かどうかは、この一次項だけからは決めない。

## 3. Nの底はcommon Y28そのものである

N normalizerの28次元sourceは、CA sourceにk1=0を課して得られる。
CA sourceの対応binary第三係数のkernel速度は零である。
moving kernel flagを戻した全150係数が、固定旗N19空間に入ることを全次数で検証する。
GL19＋残るcore9のsource augmented微分はrank28である。

Nの27独立target行はcommon Yからraw row27を除いたもの、
28番目のquadratic targetはraw row27である。
従って前稿で別に作ったNのregular R_N,28は、集合として**共通Y28の実座標環**に等しい。
ただしCA/CBの全R29との混同はしない。

同じ有限公式をN sourceに適用し、extra target row61のodd係数b_Nの一次項を求めると

    db_N=−dY16−6dY65+4dY105+dY112−4dY117

となる。quadratic branch座標row27方向の係数は零である。
Nのdiscriminant微分はそのbranch方向に非零係数−1/3を持つので、dnu_Nとdb_Nは独立である。
さらにscalar座標Y106と合わせた微分rankは3である。

## 4. projected Nの全完成像はordinary pinch point

common R29をR_N[[z61]]と見る。
N normalizerはB_N=R_N[t]/(t²−nu_N)であり、extra coordinateの像は一意に

    z61=a(Y)+b_N(Y)*t

と書ける。b_Nは上記一次項から非零である。
従ってこの一座標を加えた**実投影像**のkernelは厳密に

    (z61−a(Y))²−nu_N(Y)*b_N(Y)²

である。monic二次式で全級数を1とzへ割り、B_Nの基底1,tの独立性を使えば他のkernelはない。
nu_Nはregular parameterで非平方、b_N非零なのでこれはprimeであり、normalizerはB_Nである。

dnu_N,db_Nが独立であることから、nu_N,b_NをR_Nの二つの形式座標に選べる。
z=z61−a(Y)と置くと、**全次数で**projected Nのcone完成環は

    C[[g1,...,g26,nu,b,z]]/(z²−nu*b²)

となる。scalar Y106も独立なので、projective版はg変数25個である。
従ってこの特定の29座標射影像は滑らかなgraphではなく、ordinary pinch pointと滑らか因子の積である。
これはN自体の完成環ではない。Nの他の97 target座標がodd情報をさらに持つためである。

a(Y),b_N(Y),nu_N(Y)の全係数展開は計算していないが、既存のrank2有限source除法による
一意な全次数recipeで定まり、独立一次項だけで上のformal coordinate normal formが厳密に従う。
CAとCBの相対image kernelの和、actual Dのnilpotenceやprimary型は未解決である。

* code: `code/ambient_B18_odd_firstjets_N_projection_20261103_v63.py`
* evidence: `reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/`

全mixed source係数、二次・三次implicit補正、full odd線形行列、両discriminant微分、
Nの実source帰属とsource rank28を保存する。
