# R19境界の二枝局所環と、全normal/CM/smooth locus

2026年11月3日。R,N,B18,B17,B16は既存の全射影被約像を表す。
以下は**被約R自身と被約B18の局所環**についての結果である。
CAとCBの実scheme交叉の境界厚み、全Poisson空間の局所環はここでは計算しない。

## 結果

B17の非空稠密開上で、被約Rと被約B18の完成局所環はそれぞれ

    Rhat = C[[u1,...,u17,a,b,c,d]]/(ac,ad,bc,bd),
    B18hat = C[[u1,...,u17,s,t]]/(st).

Rhatの次元19、embedding dimension21、depth18、multiplicity2である。
normalizationは二つの正則19次元環の直積で、normalization quotientはC[[u17]]、
conductorは(a,b,c,d)である。B18hatは次元18、depth18、multiplicity2のnodeである。

さらにRにはB17の外側にもう一つnonnormalな15次元PGL5軌道C15がある。
代表点はcanonicalな

    P=∂y∧D,
    D(h)=x², D(x)=0, D(z)=hx, D(w)=h²−xz/6

である。**全Rで、点集合として**

    Sing(R) = Nonnormal(R) = NonCM(R) = B17 ∪ C15.

従ってsmooth、normal、Cohen--Macaulayの各locusは同じ開集合になる。
C15の閉包の余りはB17内に入る。
C15での完成環をcusp型と同定してはいない。

## 1. generic B17上の全normalization fiber

次を指定点とする。

    D(h)=D(x)=0, D(z)=−hx,
    D(w)=h²+2hx+3hz+4x²+5xz+6z².

intrinsic constant factorはL=C∂yで一意である。
x=0への制限でD(w)は非零なので、U∧v∧D=0を満たす定ベクトルvはまずw方向に限られ、
次にD(z)=−hxでv=0となる。ここでU=(h,x,z,w)である。
定covector Casimirsの空間はW=span(h,x)である。

既存のproper normal spectral模型で、この点の持ち上げは全てt=0上にある。
従って旗のXはW内の直線である。
Gのz²係数が非零なので、D(ell)がXで割れる線形形式の空間はK=span(h,x,z)であり、
D(z)=−hxの可除性からXはCxまたはChに限られる。
よって**全**geometric fiberはこの二点である。proper射は点の近傍で有限となり、
normal birational sourceはRのその近傍のnormalizationである。
二つのsource点はquadricのsmooth locusにある。

源のquadric tangent全11方向と全25 GL方向を、canonical150係数で計算した。
各branchのaffine rankは20、両者のunionは22である。
従ってprojectiveでは各19次元、tangent intersection17次元、全tangent21次元となる。
有限branch射の接微分が単射なので、Nakayamaにより各完成branchは正則19次元である。

B17のproper triangular模型はこの点の近傍で旗(L,W,K)を一意に戻す。
KはD²(ell)=0となる全線形形式の空間である。
その像微分のaffine rank18によりB17は指定点の近傍で滑らか17次元となる。
binary quadraticの二つの相異なる根はetale局所的に選べるので、両normalization branchは同じB17を含む。
branchのtangent intersectionが17次元であるため、両smooth branchのscheme intersectionも滑らか17次元で、
これがB17と一致する。Rはreduced excellent varietyなので完成もreducedであり、
その完成像は二つのbranchのfiber productである。
共通の17正則パラメータと各二つのnormal parameterを取ると表示された四つの積のidealを得る。

t=0のsource divisorでも同じ二枝を取り、各affine rank19、union20を直接確認した。
従ってB18の各branchは正則18次元で、同じ17次元B17に沿って横断する。
これが表示されたnode環を与える。

## 2. conductor、depthと全B17への帰結

Rの表示から正確な列

    0 → Rhat → C[[u17,a,b]] ⊕ C[[u17,c,d]] → C[[u17]] → 0

がある。最後の射は共通基底への制限の差である。
よってnormalization quotientはC[[u17]]、conductorは(a,b,c,d)、depthは18である。
multiplicityは二つの正則branchのmultiplicityの和で2となる。
B18については同じ列で各branchにnormal parameterが一つずつであり、depth18となる。

Rは既存の全域R1定理を満たす。
excellent varietyのnonnormal locusは閉で、B17は既約である。
その稠密開がnonnormalだから**全B17の各閉点がnonnormal**である。
R1にS2が加わればnormalとなるため、各点はS2を満たさず、特にCMではない。
単にfiberが連結であることからnormalityを推論してはいない。

## 3. B18\B17の全旗と、有限normalization

全B18はchain proper像なので、B17外の点は必ず

    D(h)=a x², D(x)=0, D(z)=b hx+c x², D(w)=G2(h,x,z)

でa,b≠0の表示を持つ。a=0の全像はB17、b=0の全像はB16⊂B17である。
またG2∈span(hx,x²)なら第二のCasimirを作れてB16に入る。
従ってB17外ではG2∉span(hx,x²)である。

有限座標変換

    x'=x, h'=h/a,
    z'=z/(ab)−c h/(a²b)

でD(h')=x'²、D(z')=h'x'となる。
さらにwから適当なz',h'の倍を引いてG2のhx,x²項を消せる。
全変換を物理的150係数の行列で確認した。
以降primeを省き、全点を

    D=(x²,0,hx,A h²+B hz+E xz+F z²)

と書く。(A,B,E,F)は全て零ではない。

この表示ではX=Cxが唯一のCasimir直線である。
G2 mod x≠0ならKはD(ell)がxで割れる全線形形式の空間で、span(h,x,z)に一意に定まる。
G2 mod x=0ならE≠0でD=xN_U、N_Uは四次元のregular nilpotentである。
この場合Kはその唯一の三次元invariant subspace、ker N_U³に定まる。

intrinsic Lも全ての場合に一意である。
G2 mod x≠0の場合はx=0に制限して前節と同様に示せる。
残るD=xN_Uで定ベクトルv≠0がU∧v∧N_U U=0を満たすと仮定する。
U=v+suと代入してNv∈Cvを得て、商U/Cv上のNは全直線を保つのでscalarとなる。
従ってrank(N−lambda I)≤1となるが、regular nilpotent4ではrank N=3で矛盾する。

したがって全B18\B17においてproper spectral模型のgeometric fiberは一つの旗に尽きる。
sourceはd=1のquadric chartにあり、smooth19次元である。
propernessで有限化した後、このsourceはtargetの**有限normalization**である。

## 4. 全ramificationの正確な線形計算

quadric tangentの全11列とGLの全25列による150×36行列をJとする。
19個の列と行を取るminorは**全A,B,E,Fで定数48**である。
その列を消去すると、残る行列の非零entryは一つの列にのみあり、正確に

    B, A/6+E, F, −B, −A/6−E, −2F

である。全entryと全minorを証明書に保存した。
従ってsource differentialのaffine rankは

    20  iff (B,F,A/6+E) != (0,0,0),
    19  iff B=F=A/6+E=0.

前者はfinite normalizationのfiberがreduced length1なので、Nakayamaで射が局所同型となる。
従ってRはそこで滑らかである。
後者はB17外なのでA≠0であり、wをAで割ると全て指定されたC15の代表点になる。
そのGL-affine orbit rankは16、従ってprojective軌道次元は15である。
この場合normalization differentialは単射でないので射は局所同型ではない。
もしtargetがnormalなら有限birational射は同型となるはずであるため、target Rはnonnormalである。
全域R1によりここもnonCMである。

既存定理でR\Nは全てsmoothであり、残る全境界はB18である。
上の全chart exhaustion、全B17のnonnormality、C15の一軌道分類を合わせ、
記載した三つのlocusの正確な点集合等式を得る。
閉性からC15の閉包の追加点はB17内に限られる。

## 証明書と範囲

* `code/ambient_R19_boundary_local_models_20261103_v63.py`
* `reports/AMBIENT_R19_BOUNDARY_LOCAL_MODELS_20261103_v63/report.json`

全QQ計算1.95秒。二つのR branch・B18 branch・union・B17の全150行列とminor、
全係数で定数48のminorおよび消去残差、全chainの有限座標変換を保存した。
C15の完成環の生成関係やB16のquadric-vertex上の局所環は未計算である。
実CA/CB交叉のboundary primary structureをこの被約像から推論しない。
