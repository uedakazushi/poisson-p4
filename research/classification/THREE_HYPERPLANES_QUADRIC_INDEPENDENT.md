# 独立な三超平面と既約二次曲面：全配置、非零対数括弧の被覆、新23次元成分

2026-11-07。課題 E_reduced_nonnormal。
実退化因子を D=H0+H1+H2+Q とし、三つの超平面は相異なり一次独立、Qは既約二次曲面とする。
同一超平面の重複度三を意味する3H+Qとは異なる課題である。

本稿では全配置を有限個に帰着し、非零対数括弧の全支持を既知MPS40成分の閉和で覆う。
対数括弧が零の枝は四つの配置と各全六磁気係数へ帰着する。
三辺全てに接する配置の全像は新しい23次元既約成分であり、残る三配置の全六磁気源も実変形によりMPSへ覆う。
従ってこの実退化因子型の独立三超平面の場合は全支持が被覆される。

## 1. 全二次模型と二つの枝

三超平面をh=0,x=0,y=0とする。実因子hを使うh-Casimir gaugeにより、
h=1上の二次bivector P はstrict Jacobiを満たす。
x,yも実因子なのでnormalであり、全括弧は

    {x,y}=a xy,       {x,z}=x A,       {x,w}=x B,
    {y,z}=y C,        {y,w}=y D,       {z,w}=F,

となる。aは定数、A,B,C,Dは任意のaffine一次式、Fは二次式である。
q=Q|_(h=1)とすると、全Pfaffian条件は

    aF−(AD−BC)=kq,               k≠0.                 (1)

以下の議論はqの次数落ちやQの階数三・四・五を落とさない。
a=0とa≠0は三超平面の置換で不変である。
例えばchart hからchart xへ移ると対数座標は(log x,log y)から(−log x,log y−log x)へ移り、aの符号だけが変わる。

## 2. a≠0：可換な二つの線形場による全表示

V=A∂z+B∂w、W=C∂z+D∂w、X=x∂x、Y=y∂yと置き、

    D1=X−W/a,        D2=Y+V/a,        b=k/a

と定義する。全Jacobiと(1)は正確に

    [D1,D2]=0,
    [D1,q∂z∧∂w]=[D2,q∂z∧∂w]=0,
    P=aD1∧D2+bq∂z∧∂w,             ab≠0              (2)

と同値である。
逆方向も成立し、二次元fiberのbivectorは自動的にPoissonなので、(2)は全Jacobiを保証する。
D1,D2はhを固定する斉次線形場に延びる。
第二式は

    Di(q)=tr_vertical(Di) q

という相似保存条件である。

### 相対配置の有限化

二つの作用のjoint characterを使い、h,x,yのcharacterをそれぞれ0,e1,e2と書く。
e1,e2は独立である。qのcharacterをτ、quadricのradicalをRとすると、Rは両作用で不変であり、

    tr_R(D1,D2)=e1+e2+(1−rank(Q)/2)τ.                 (3)

これは全五次元traceがe1+e2+τであり、非退化商のtraceがrank(Q)τ/2であることから従う。
Ann(R)上の逆二次形式に対し、固有covectorのcharacter和がτと異なればそのpairingは零である。
ここでは両作用の半単純性を仮定しない。

rank(Q)=5なら(3)はτ=2(e1+e2)/3を与える。
h,x,yの全pairのcharacter和はこの値と異なるため、三次元の全等方部分空間ができる。
五次元非退化二次形式のWitt indexは二なので矛盾する。従ってrank五は空である。

rank四ではRは一次元であり、h,x,yのうちR上で非零なものは高々一つである。
三超平面が全て頂点を含む場合、Ann(R)内の三次元Gram行列のrankは少なくとも二。
character比較から、その非零部分はh-x、h-y、x-yのいずれかの一pairとなる。
三超平面を置換し、z,wの線形変更とh,x,yのfiberへのshearを用いるとq=hz+xyとなる。
ちょうど二超平面が頂点を含む場合、残る超平面をhとし、頂点をh-axisへ送る。
qはx,y,z,w上で非退化、τ=e1+e2、x,yは共に等方covectorである。
その相互pairingが非零か零かにより、q=xy+zwまたはq=xz+ywとなる。

rank三ではRは二次元であり、h,x,yのRへの非零restrictionは高々二本。
三超平面全てがvertex lineを含む場合、qはh,x,yだけの式となる。
六つの二次単項式のcharacterは相異なるため、相似保存されるqは単一単項式となり、既約性に反する。
ちょうど二超平面がvertex lineを含む場合、それらをx,yとし、R=span(h,w)へ送る。
三次元非退化商内でx,yのGram rankは少なくとも一。
τ=2e1,e1+e2,2e2のいずれかとなり、q=x²+yzまたはq=xy+z²を得る。
ちょうど一超平面がvertex lineを含む場合、それをhとしR=span(x,y)へ送る。
(3)からτ=0。hの接触・非接触に従いq=hz+w²またはq=h²+zwとなる。

従ってa≠0の全点は、三超平面の置換を含め、次の七配置で覆われる。

| rank Q | q | vertexの配置 |
|---|---|---|
| 4 | xy+zw | 二超平面だけが頂点を含む |
| 4 | xz+yw | 同上、二つの接hyperplaneのdual pairが零 |
| 4 | hz+xy | 三超平面全てが頂点を含む |
| 3 | h²+zw | hだけがvertex lineを含み、hは非接触 |
| 3 | hz+w² | hだけがvertex lineを含み、hは接触 |
| 3 | xy+z² | x,yがvertex lineを含み、両方接触 |
| 3 | x²+yz | x,yがvertex lineを含み、一方だけ接触 |

### 二つの全源の直接SNC極限

q=xy+zwでは、相似保存の全係数比較から

    D1=x∂x+α z∂z+(1−α)w∂w,
    D2=y∂y+β z∂z+(1−β)w∂w

が全解である。h方向のfiber translationも非対角項も残らない。
重み(h,x,y,z,w)=(0,1,1,0,0)で(2)の全項は非負重みとなり、
正重みのbxy∂z∧∂wだけが消える。
零重み部分はlogcanonicalで、その退化式はab hxyzw、従って被約SNCである。

q=h²+zwでも全解は

    D1=x∂x+α(z∂z−w∂w),
    D2=y∂y+β(z∂z−w∂w).

今度は重み(1,0,0,0,0)でbh²∂z∧∂wだけが消え、同じ非零SNC退化式を得る。
この二つは全α,βおよび全ab≠0についての極限である。

### 残る五配置の全nilpotent源

相似保存と可換性を全係数で比較し、qを保つfiber shearを施すと、各配置の全源は

    D1=D1⁰+pN,       D2=D2⁰+qN                         (4)

となる。ここで(4)のscalar qは表のquadricと別の自由係数であり、p,qは全複素数を動く。

| quadric | D1⁰ | D2⁰ | N |
|---|---|---|---|
| xz+yw | x∂x+w∂w | y∂y+z∂z | y∂z−x∂w |
| hz+xy | x∂x+z∂z | y∂y+z∂z | h∂w |
| xy+z² | x∂x+(z∂z+w∂w)/2 | y∂y+(z∂z+w∂w)/2 | z∂w |
| x²+yz | x∂x+2z∂z | y∂y−z∂z+w∂w | y∂w |
| hz+w² | x∂x | y∂y | h∂w−2w∂z |

全源性を確かめるため、shear前の余剰項も記す。
hz+xyではD1の∂w係数はph+rx+sz、D2ではqh+ty+szであり、
wをw−rx−ty−szへ替えれば(4)となる。
xy+z²では共通のh係数、互いに符号反対のx係数とy係数だけが非共鳴項であり、
これら三項をwの一次shearで除くと、z∂wの二つの独立係数が残る。
x²+yzではD1の追加係数はrx+py+sz、D2ではuh−rx+qy−sz。
wのh,x,z shearでr,s,uを除くとy∂wの二係数が残る。
xz+ywとhz+w²では、表示したN以外の係数は相似保存そのものから零になる。
従って共鳴係数を零に置いて稠密部分だけを扱っているのではない。

### 全源を保つ実Poisson変形

各行に対し、次のqε,Nεを用いる。

| q0 | qε | Nε | ε≠0の既被覆配置 |
|---|---|---|---|
| xz+yw | xz+yw+εzw | (y+εz)∂z−(x+εw)∂w | rank4、xy+zw型 |
| hz+xy | z(h+εw)+xy | (h+εw)∂w−εz∂z | rank4、xy+zw型 |
| xy+z² | xy+z²+εw² | z∂w−εw∂z | rank4、xy+zw型 |
| x²+yz | x²+yz+εzw | (y+εw)∂w−εz∂z | rank3、h²+zw型（超平面を置換） |
| hz+w² | hz+w²+εz² | (h+2εz)∂w−2w∂z | rank3、h²+zw型 |

各行で、Nε(qε)=0、div_vertical Nε=0、[D1⁰,Nε]=[D2⁰,Nε]=0が恒等的に成立する。
D1⁰,D2⁰もqε∂z∧∂wを保存する。
従って任意のp,q,a,b（ab≠0）について

    Pε=a(D1⁰+pNε)∧(D2⁰+qNε)+bqε∂z∧∂w             (5)

は全εにわたる真正の二次Poisson族である。
ε=0で任意の正規化済み元へ戻り、必要なら先ほどの固定fiber shearを逆に戻せる。
ε≠0では上表の二つの直接SNC極限を持つ全源へ入る。

MPS Theorem1.4と既存40成分の同定を外部入力として使う。
直接SNC極限を持つ元を含む各全Poisson成分は、PGLの連結性によりそのSNC極限も含み、MPS成分である。
既知40成分の有限和は閉じているため、実族(5)のε=0もその和に入る。
これでa≠0の七配置の**全支持**が既知MPS40の閉和に含まれる。
特定の一つのMPS番号や、元を通る全成分の一意性は要求していない。

## 3. a=0：rank三のbase conicと四配置

Jacobiから[V,W]=0となる。fiber変数をZ=(z,w)とし、
V=AZ+u(h,x,y)、W=BZ+v(h,x,y)と書けば、

    [A,B]=0,       Av=Bu.

span(A,B)に可逆行列があれば、V,Wを一定の可逆結合で替えてAを可逆とできる。
するとW=BA⁻¹Vであり、det(V,W)はVの二成分だけの二次形式となって複素数体上で因数分解する。
span(A,B)が全てsingularで片方が非零なら、2×2可換行列の中心化を用いてB=cAとなる。
Aが非零固有値を持つrank一型かnilpotent rank一型かに従い、Av=Buは
det(V,W)を二つの一次式の積にする。
(1)とQの既約性に矛盾するので、必ずA=B=0である。

従ってV,Wはbase(h,x,y)だけに依るtranslationで、Q=q(h,x,y)はrank三。
残る二つのJacobi式ではVF,WFがbaseだけの函数となる。
det(V,W)≠0の函数体上で逆を取るとFのfiber二次部分が零となり、Fはfiber-affineである。
この結論は元の多項式環へ戻る恒等式であり、Q=0上の点を除外していない。

L={h=x=y=0}を除いた射影空間はTot(O_P²(1)⊕O_P²(1))である。
fiber-linear部分はrank二Lie algebroid E=O_P²(−1)²を与え、そのanchorの退化因子は

    B={hxyq(h,x,y)=0}⊂P².

JacobiによりanchorはBに対数的に接する。
被約plane curveの対数接束は局所自由であり、
det E=O(−2)=det T_P²(−log B)なので、anchorはこの二つの束の同型となる。
従ってBはexponents(2,2)のfree plane quinticであり、全Tjurina数は12である。

三角形の頂点のうちconic上にある数をm、conicに接する辺の数をtとする。
頂点での接触もtに一回数える。
可能な局所型はA1,A3,D4,D6であり、全pairwise intersection multiplicityが9なので

    τ(B)=9+m+t.

従ってm+t=3。
各conicは滑らかであるから、三角形を保つ対角変更と辺置換の下で、可能な配置は次の四つに尽くされる。

| (m,t) | qの正常形 | 接触の配置 |
|---|---|---|
| (0,3) | h²+x²+y²−2hx−2hy−2xy | 三辺全てに接する |
| (1,2) | (x+y)²+hy | 既存頂点で一つのincident辺に接し、opposite辺にも接する |
| (2,1) | hx+xy+y² | 二頂点を結ぶ辺以外の一辺で既存頂点に接する |
| (3,0) | hx+hy+xy | 三頂点を通り、接する辺はない |

完全性はconicの六係数の比較でも確認できる。
三頂点を通る場合は三つの対角係数が零、二頂点の場合は一つの交差係数も零、
一頂点の場合はopposite辺の判別式が零となる。
頂点を通らず三辺に接する場合は三つの辺判別式が零となり、
rank一の平方を除く符号だけが残る。各行の非零係数は対角変更で上表へ規格化できる。

固定した各配置とanchorの基底を選ぶと、Lie bracketがFのfiber-linear部分を一意に決める。
残る自由度は

    F0(h,x,y)∂z∧∂w,       F0∈Sym²⟨h,x,y⟩,

という**全六係数**である。
これは任意のbase函数のHamiltonian場がverticalであるため、元のPと常にcompatibleである。
六係数を標本に固定したり、無断で除いたりしない。
全六係数を保持した有限計算と成分性については第5節で扱う。

## 4. 単なる軌道極限では残余を閉じられない理由

三つの独立超平面の共通交線Lに対し、rank(Q|L)≤1という条件は、
L∈Gr(2,5)も保持した三超平面・二次式・Grassmannianのproperな積空間上で、
Hi|L=0およびrank(Q|L)≤1という閉incidenceで表される。
これにより三超平面がdependentになる境界も落とさない。
SNCの五独立超平面への極限では、選ばれた三超平面はなお独立で、
極限のLはその共通交線と一致する。
残る二つの一次式のLへのrestrictionは独立なので、積二次式のrankは二になる。
従ってrank(Q|L)≤1の元には、PGL軌道だけによるSNC極限は存在しない。
第2節の残る五配置では、この障害を実Poisson族(5)によって越えた。
第3節の四配置でも、同じ障害を無視してSNC所属を推論しない。

## 5. 全六磁気係数、四つの23次元像、および一つの新成分

affine h=1で固定したanchor基底をU,Vとし、[U,V]=cU+dVと書く。
c,dはaffine一次式であり、全源は

    P_R=U∧∂z+V∧∂w+(−cz−dw+R(x,y))∂z∧∂w,
    R∈C[x,y]≤2.                                             (6)

全六係数でJacobiは恒等的に零となり、Pfaffianは−xyqである。
二次syzygyの全18係数を解いて得る二次元核からU,Vを取るため、
anchorや磁気係数の一部を標本で置換していない。
別の基底はfiberの定数GL2変更で吸収され、(6)は固定配置の全源を与える。

各源のPGL飽和像の次元は**23**である。
上界は次の全パラメータの対称性から従う。
PGL5の七次元部分群

    z′=ρz+φ(h,x,y),   w′=ρw+ψ(h,x,y),
    ρ≠0,   φ,ψ∈⟨h,x,y⟩

は、射影スカラーρを除いた後、

    R′=ρR+(V+c)φ+(d−U)ψ                                  (7)

により(6)を保存する。
PGL5×A6からPoisson空間への写像の各fiberは、この部分群の自由な七次元軌道を含む。
従って像次元は24+6−7=23以下。
下表の実source微分階数48は、Euler verticalの24方向とスカラー一方向を除いて
23という下界を与えるため、全四像について上界と一致する。

全150係数・全350射影Jacobi式を使った固定証拠点の結果は次の通りである。
実sourceにはPGL、スカラー、全六磁気微分を含めた。

| 配置 | Dの階数 | 実source階数 | 源像次元 | 射影接次元 | 真正法線数 |
|---|---:|---:|---:|---:|---:|
| m=0,t=3 | 102 | 48 | 23 | 23 | 0 |
| m=1,t=2 | 101 | 48 | 23 | 24 | 1 |
| m=2,t=1 | 100 | 48 | 23 | 25 | 2 |
| m=3,t=0 | 99 | 48 | 23 | 26 | 3 |

主証明書は
`reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_20261107.json`、
独立証明書は
`reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_REVIEW_20261107.json`。
独立実装は全syzygy核、普遍Jacobi、Pfaffian、四つの全行列、実sourceと法線による
全接核をQQで再構成している。

### 5.1 m=0,t=3の新23次元成分

この配置では、より簡単な同値な基底を

    s=x+y−h,
    U=(xs,2xy),   V=(2xy,ys)

と取れる。h=1で[U,V]=−yU+xVなので、

    P_R=xs∂x∧∂z+2xy∂y∧∂z+2xy∂x∧∂w+ys∂y∧∂w
        +(yz−xw+R(h,x,y))∂z∧∂w                         (8)

が全六磁気係数の正常形である。
退化式は

    D=−hxy{h²+x²+y²−2hx−2hy−2xy}.

PGL5×A6は既約で、その像閉包Zも既約。
上表の証拠点では全空間の射影接次元23とdim Z=23が一致する。
従ってその点は全Poissonスキームの滑らかな点であり、Zは23次元の全空間既約成分である。
この論法は残る三配置の成分性を含意しない。

Zは既知の成分と相異なる。
まず五次式の本質変数数が高々三である条件は、五偏微分の係数行列の階数条件として閉である。
Z全体でこの条件が成り立つが、被約SNC五次因子では成り立たない。
実際、三変数五次式の二次以上の既約因子は共通の頂点P¹で特異となり、
全て一次因子なら五超平面がそのP¹を共有してSNCにならない。
従ってSNC点を持つMPS40成分のいずれとも一致しない。
この閉包全体での比較は
`ADDITIONAL_RANK4_LATE_COMPONENT_AUDIT.md`第3節と同じである。

一般rank四なので六つのrank二成分と異なり、次元により22、24、28次元の追加成分とも異なる。
既知23次元追加成分のうち、H_QUARTIC、M、A0、DPM、ZERO_LOG_5Hは
一般退化式の因数分解型が異なる。
残るCMINUS_MAGの一般因子は三超平面と既約二次式であるが、三超平面はdependentである。
dependent三超平面を持つ条件はproperな因子incidenceの像として閉で、
(8)の独立なh,x,yを持つ点はそこに入らない。
従ってZはこれらからも異なる新成分である。

### 5.2 m=1,2,3の全磁気源をMPSへ入れる実変形

ここでは接空間や有限jetから実変形の存在を推測せず、全パラメータで定義されたPoisson族を与える。
以下のt,α,β,γは独立なパラメータである。まず座標(h,x,y,Z,W)において、

    D1=x∂x+p Z∂Z+ε1 W∂W,
    D2=y∂y+r Z∂Z−r W∂W

という可換な対角場を取る。表の二次式QはD1,D2についてそれぞれ
vertical traceと同じcharacterを持つので、[Di,Q∂Z∧∂W]=0である。
従って

    Π_t=t D1∧D2+b Q∂Z∧∂W                              (9)

は全パラメータのstrict Poisson構造となる。

| 配置 | ε1 | p | r | Q | b |
|---|---|---|---|---|---|
| m=1,tangent数2 | 1−p | tα | 1/2+tβ | ZW/2−hx | −t(1+tγ) |
| m=2,tangent数1 | −p | −1+tα | 2+tβ | ZW+h² | −2t(1+tγ) |
| m=3,tangent数0 | −p | 1+tα | −1+tβ | ZW−h² | t(1+tγ) |

この表のtは変形パラメータであり、従来の接線本数tとは区別する。
以降は配置をmだけで表示する。
次の座標変更で(9)を(h,x,y,z,w)へ引き戻す。

| m | Z | W |
|---|---|---|
| 1 | t z+2(x+y) | t(w−z/2)+h+x+y |
| 2 | t z+h+y | t w+x+y−h |
| 3 | t z+h+x | t w+h+y |

各座標変更の行列式はt²。
引き戻した式は見かけ上t⁻¹を持つが、その分子の定数項は全パラメータで零となる。
従って実際のΠ_tはtを含む**多項式**の二次bivectorであり、t=0にも正則に延長する。
証明書では全十個のstrict Jacobi式をこの多項式から直接再検査している。

t(1+tγ)≠0では、対数括弧の係数はa=t≠0。
Qの階数はm=1で四、m=2,3で三であり、第2節で全支持をMPSへ被覆済みの配置になる。
従ってこの開集合上の全点はMPS40成分の和に入る。

t=0のmixed anchorとbracket係数は次の通りである。
これらは第3節の各free配置の二つの基底であり、定数GL2変更により元の基底と同値となる。

| m | U | V | c,d |
|---|---|---|---|
| 1 | (x(x−y),2xy) | (−x(h+4y)/2,−y(h−x+y)) | h+y,x |
| 2 | (x(2h+y),y(h+y)) | (x(2h−2x−3y),y(h−y)) | y−2h,2h+y |
| 3 | (−x(h+x),−hy) | (hx,y(h+y)) | −h,−h |

すなわち各極限は(6)と同じ形であり、磁気項は

    R=αRα+βRβ+γRγ

となる。その三つの係数と、fiber translationから得るcoboundaryの像は次である。

| m | Rα | Rβ | Rγ | translationの像 |
|---|---|---|---|---|
| 1 | −2hy−4xy−4y² | 2x²−2y²−2hy | −(x+y)²−hy | ⟨h²+hy,hx,xy⟩ |
| 2 | −xy−2y² | hx+xy | −2(hx+xy+y²) | ⟨h²,hy,x²+xy⟩ |
| 3 | −hy−xy | hx+xy | hx+hy+xy | ⟨h²,x²,y²⟩ |

各行で三つのR係数と右端の三つの基底はSym²⟨h,x,y⟩全体を張る。
係数順を(h²,hx,hy,x²,xy,y²)、列順を表の三つのR、右端の三基底とすると、
六次正方行列式はそれぞれ8,2,−1で非零である。
従って任意の六係数磁気項が、(9)のt=0極限と射影的fiber translationから得られる。

MPSの有限個の閉成分の和は閉であり、PGL5で不変である。
以上により**m=1,2,3の全六磁気源**がMPS40成分の和に含まれる。
α,β,γの特殊値や磁気項零の場合を除外していない。
t(1+tγ)=0の補集合もこの閉性で保持される。

全多項式族、正確な極限、全Pfaffian、全coboundaryおよび磁気係数行列を含む自足証明書は

- `reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_M1_MPS_20261107.json`
- `reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_M23_MPS_20261107.json`

に保存した。主実装は対応する`code/classification_three_hyperplanes_quadric_m1_mps.py`と
`code/classification_three_hyperplanes_quadric_m23_mps.py`であり、各reportにも実装全文を固定している。
別担当は三つのpole cancellation、極限anchor、三つの磁気係数、全coboundary像および
六係数到達を手計算で独立に確認した。
さらに`reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_MPS_REVIEW_20261107.json`は、
別のFraction疎多項式実装で全十Jacobi、実Pfaffian、全極限、
座標変換後のcanonical恒等式TΠTᵗ、二次式階数、旧全源との定数GL2基底変更、
全coboundaryおよび三つの六次行列式を独立に再構成している。

この被覆証明は、その前に保存した二次・三次法線診断の成分性推測には依存しない。

## 6. 独立三超平面の全支持

第2節はa≠0の全配置をMPSへ覆い、第3節はa=0の全配置を四つに尽くした。
第5節はm=0の全源を新23次元成分とし、m=1,2,3の全源をMPSへ覆った。
従って、相異なる独立三超平面と既約二次曲面を実退化因子に持つ**全Poisson構造**は、
MPS40成分と第5.1節の新23次元成分の和に含まれる。
三超平面の従属配置は別ノートの定理と合わせて扱う。

## 独立確認と射程

rank五の排除は、三つの既知固有covectorだけを用い、半単純性を仮定しないことを独立に確認した。
a=0の四配置、頂点接触を含むτ=9+m+t、各正常形の完全性は別担当が読解確認した。
第2節の全線形場の係数比較、全nilpotent係数を保つ(5)、二つのSNC極限は手計算で独立に照合した。
MPSの外部定理自体は既存の依存である。

三超平面が二次元のpencilを張る場合は別ノートで扱う。
全Poisson空間の分類、および同一超平面を重ねた3H+Qへの拡張は主張していない。
