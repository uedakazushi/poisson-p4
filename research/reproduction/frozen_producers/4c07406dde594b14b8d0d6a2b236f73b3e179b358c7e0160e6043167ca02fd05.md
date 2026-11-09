# 独立な三超平面と既約二次曲面：全配置と非零対数括弧の被覆

2026-11-07。課題 E_reduced_nonnormal。
実退化因子を D=H0+H1+H2+Q とし、三つの超平面は相異なり一次独立、Qは既約二次曲面とする。
同一超平面の重複度三を意味する3H+Qとは異なる課題である。

本稿では全配置を有限個に帰着し、非零対数括弧の全支持を既知MPS40成分の閉和で覆う。
対数括弧が零の枝は四つの配置と各全六磁気係数へ帰着するが、その既知帰属は未証明である。

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
この四配置の全六係数源の既知成分への帰属は、次の未解決義務として保持する。

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

## 独立確認と射程

rank五の排除は、三つの既知固有covectorだけを用い、半単純性を仮定しないことを独立に確認した。
a=0の四配置、頂点接触を含むτ=9+m+t、各正常形の完全性は別担当が読解確認した。
第2節の全線形場の係数比較、全nilpotent係数を保つ(5)、二つのSNC極限は手計算で独立に照合した。
MPSの外部定理自体は既存の依存である。

三超平面が二次元のpencilを張る場合は別ノートで扱う。
全Poisson空間の分類、残るa=0四源の包含、同一超平面を重ねた3H+Qへの拡張は主張していない。
