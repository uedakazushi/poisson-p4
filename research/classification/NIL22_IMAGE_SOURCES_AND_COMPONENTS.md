# nil22 image内：全十一源と二つの23次元成分

2026-11-07。`E_multiple_hyperplane`。重複超平面を無限遠とし、最高二次部分の
Pfaffian が零、canonical affine modular 場が Jordan (2,2)、定数部がimage内にある
**全源**を扱う。全源は十一個の明示的既約多項式源へ帰着する。
そのうち `triple_A`,`triple_B` の全六磁気源は相異なる23次元全空間成分を与える。
さらに `quadratic_B`,`triple_C` の全九変数源、`double_D` の全八変数源は、
実Poisson曲線で既済half/third/CN31の25成分へ入る。
残六源の成分帰属はここでは未解決であり、有限表示だけを被覆完成とは呼ばない。
先行文献上の新規性は未確認。統合担当による計数は別作者reviewのPASS後に限る。

## 1. 全源の必要十分な二変数行列表現

座標を `(q1,q2,x,y)`、`p=(x,y)`、モジュラー場を
`N=x∂q1+y∂q2` とする。`MODULAR_NIL_RANK2_FIRSTJET_BOUND.md` §2 の
全90係数reductionにより、全点は

    P=F∂q1∧∂q2+ΣMij(p)∂qi∧∂pj,
    F=F0(p)+(y+div M2)q1+(-x-div M1)q2

である。`F0=m0+m1x+m2y+m3x²+m4xy+m5y²` は任意で、六係数を全て保持する。
`D=det M≠0`、`degD≤3`、`D=D0+D2+D3`、`S=Hess(D2)` と書くと

    H=S+Hess(D3)/2,       Hp=∇D,
    M=-adj(H)+p⊗V(p),     V=V0+V1p,
    D0=det S,
    D=det H−V·∇D.                                      (1)

ここまでは同ノートの必要条件である。本稿ではこれが十分であることも使う。
実際、(1)より `M∇D=-Dp`。Mの二行をベクトル場U,Vとすると、D≠0の開で

    [U,V]=(div V−V(logD))U+(U(logD)−div U)V
         =(div V+y)U+(−x−div U)V.

従って `Fq1=divV+y`,`Fq2=−x−divU` が二つの非自明Jacobi式をちょうど満たす。
残る二式はPpp=0とMのq独立性で零。任意F0も許される。
全て多項式恒等式なのでD=0上にも延びる。
従って(1)はrank4開に対する必要十分条件である。

同時変換 `(q',p')=(Aq,Ap)` はNを保ち、

    M'=AMAᵀ,    D'(p')=det(A)²D(A⁻¹p'),
    S'=det(A)²A⁻ᵀSA⁻¹,   V'=A V(A⁻¹p').              (2)

F0は可逆線形変換で任意二次式全体へ移る。以下の正規化で磁気係数を捨てない。
複素数上、非零binary cubicは三単根・二重根・三重根の三型に分かれ、(2)で
`xy(x−y),x²y,x³` へ正規化できる。D3=0も別に保持する。

## 2. 三単根：二つの全六磁気源

D3=xy(x−y)では最高次数式よりV1=−I/3。
`S=[[u,v],[v,w]]` として全低次係数を比較すると

    V0=(1+w/6,−1−u/6),  v+u+w+9=0,
    f1=u²+2uw+15u−6w−54=0,
    f2=w²+2uw−6u+15w−54=0.

`f1−f2=(u−w)(u+w+21)`。
u=wでは `f1=3(u+6)(u−3)`、u+w=−21では `f1=−(u−3)(u+24)` なので、
全Hessianは

    [[3,−15],[−15,3]], [[−6,3],[3,−6]],
    [[3,12],[12,−24]], [[−24,12],[12,3]]

の四個。Hessianの可逆性は仮定せず、singular補集合もこれで空と分かる。
一番目から四番目へのA=`[[1,0],[1,−1]]`、四番目から三番目へのA=`−swap`
で全Mと全六磁気源が一致する。従って二源を残す。

|源|D=detM|
|---|---|
|stable_line_conic|(x−y−12)(2xy+3x−3y+36)/2|
|stable_three_lines|(x+3)(y−3)(x−y−3)|

完全計算は `CLASSIFICATION_NIL22_IMAGE_STABLE_CUBIC_20261107.json`。
両源の全D階数99、実F階数48、真正法線数3。
最初の真正二次支持はともに正次元であり、成分性・既知族帰属はまだ従わない。
この失敗から高次jetの自動反復は行わない。

## 3. 二重根：四つの源と全補集合

D3=x²yでは `V1=diag(a,−1−2a)`。同じSを用いると

    V0=(a v/2,−1−(a+1/2)u),
    (4a+1)w=0,
    v[2−(a+1)u]=0,
    av²−[4+(2a+1)u]w=0.                               (3)

w≠0ならa=−1/4。
v≠0ではu=8/3,w=−3v²/64、v=0ではu=−8。
w=0,v≠0ではa=0,u=2。
w=v=0ではa,uが任意。
二重根を保つA=diag(s,1)によりv→sv,w→s²wなので、各非零vまたはwを1へ正規化する。
零補集合は全て(3)の別枝に残っている。

|源|S|V0|V1|追加自由変数|D|
|---|---|---|---|---|---|
|double_A|[[8/3,1],[1,−3/64]]|(−1/8,−5/3)|diag(−1/4,−1/2)|なし|(384x²y+512x²+384xy−9y²−432)/384|
|double_B|diag(−8,1)|(0,1)|同上|なし|(y−4)(2x²+y+4)/2|
|double_C|[[2,1],[1,0]]|(0,−2)|diag(0,−1)|なし|(x+1)(xy+x−1)|
|double_D|diag(u,0)|(0,−1−(a+1/2)u)|diag(a,−1−2a)|a,u|x²(y+u/2)|

## 4. 三重根：三つの源と全補集合

D3=x³では `V1=[[-1/3,0],[c,d]]`,`V0=(m,n)`、全式は

    m=−u/18−cv/3,
    (d+2/3)v+cw=0,   (d+1/2)w=0,
    um+vn=3w,        vm+wn=0.                          (4)

w≠0ではd=−1/2,c=−v/(6w)。Δ=uw−v²とすると `Δ²+54w³=0`、特にΔ≠0。
D3を保つ `A=[[r²,0],[t,r]]` では
`w'=r⁴w`,`v'=r³v−r²tw`,`u'=r²u−2rtv+t²w`。
t=rv/w、r²=−18w/Δを選ぶとS'=diag(−18,−6)、V0'=(1,0)、
V1'=diag(−1/3,−1/2)となる。
w=0,v≠0ではm=n=0,d=−2/3,c=−u/(6v)。同じ群でr³v=1、t=ru/(2v)を取り、
S'=[[0,1],[1,0]]、V0'=0、V1'=diag(−1/3,−2/3)とする。
w=v=0では(4)からu²=0。被約支持上u=m=0、n,c,dは任意。

|源|S|V0|V1|追加自由変数|D|
|---|---|---|---|---|---|
|triple_A|diag(−18,−6)|(1,0)|diag(−1/3,−1/2)|なし|x³−9x²−3y²+108|
|triple_B|[[0,1],[1,0]]|0|diag(−1/3,−2/3)|なし|x³+xy−1|
|triple_C|0|(0,n)|[[-1/3,0],[c,d]]|n,c,d|x³|

## 5. D3=0：二つの源

D=detS+pᵀSp/2、`SV0=0`,`V1(D2)=−D2`。
S非退化なら複素GL2でS=[[0,1],[1,0]]へ正規化し、V0=0、
V1=diag(λ−1/2,−λ−1/2)。SのS-skew endomorphismが一次元であるためこれが全て。
rankS=1ならS=diag(1,0)、V0=(0,n)、V1=[[-1/2,0],[c,d]]。
S=0ならD=0でありrank4点を持たない。

|源|追加自由変数|D|
|---|---|---|
|quadratic_A|λ|xy−1|
|quadratic_B|n,c,d|x²/2|

§3–5の九源は `CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json`。
全源について(1)、全Jacobi、divP=N、最高Pf(Q)=0を自由パラメータの多項式恒等式として検証した。
§2と合わせ十一源は元の全rank4被約支持を尽くす。各源自体は既約affineパラメータ空間であるが、
異なる源の閉包が相異なるという主張は、以下の二源以外には含めない。

## 6. triple_Aとtriple_Bの真正切片支持

固定磁気係数はm=(1,2,3,5,7,11)。各源の分母を払う整数tensor倍で全150係数を固定し、
全350のprojective Jacobiを微分する。Euler24、PGL、scalarと全六磁気微分を含む実Fを用いる。

|源|rankD|rankF|dim kerD|真正法線数|射影接次元|
|---|---:|---:|---:|---:|---:|
|triple_A|101|48|49|1|24|
|triple_B|100|48|50|2|25|

DF=0、DN=0、rank[F,N]=dimkerDをQQ上で確認。
Dの独立な101または100本の標準座標列に真正法線を足した102次元Lについて
rank[F,L]=150を確認した。従ってLは実際の48本の線形cutに対する真の切片である。
全D150列を消すfunctionalを構成し、normal二次Jacobiの全pairingを取ると、
適切な有理正規化後の初期式は

    triple_A:  a²,
    triple_B:  q1=a²−(3/124)b²,  q2=ab+(5/31)b².

B=5/31,A=3/124とすれば

    B²−A=7/3844≠0,
    b*q1−(a−B*b)*q2=(B²−A)b³.

従って双方の根基はnormal原点。
独立D式の形式陰関数による消去後、functionalが全D像を消すため上記は残idealの真正初期式となる。
その最大根基性から真の102次元切片の被約支持は原点であり、切片Krull次元は零。
48本のcutの次元不等式により全150係数空間の局所次元は48以下。
実F48からの下界と一致するので、各全六磁気源のPGL飽和閉包は全空間成分であり、
Euler24と射影scalar1を引くと23次元となる。
滑らかさ、局所成分一意性、正確なnil指数、形式積分解は主張しない。

## 7. 二成分の相違と既存94成分の監査

全源の退化五次式は非零scalarを除いてそれぞれ

    DA=h²(x³−9hx²−3hy²+108h³),
    DB=h²(x³+hxy−h³).

六磁気係数には依存しない。両三次式は既約でnodeを持つ。
DAのaffine三次式は `3y²=(x−6)²(x+3)` と書け、右辺はC(x)の平方でないので
C上で既約。DBはyに関して一次で、係数hxと定数x³−h³がC[h,x]で互いに素なので
Gauss補題でC上既約である。
DAのnodeは[h:x:y]=[1:6:0]、唯一の二重超平面h=0はnodeを通らず、
三次曲線の滑らかなflex [0:0:1]で三重接触する。
DBのnodeは[0:0:1]、その局所接錐はhx、h=0はnodeでの一方の接線。
唯一の二重超平面と唯一の既約三次因子は内在的なので、このincidenceはPGL不変である。
同じ23次元成分なら両PGL飽和constructible稠密像は共通の非空開を含み矛盾する。
従って相異なる二成分である。

全閉包はproperな掛け算写像の閉条件「D=H²F3」に入る。
既存23次元成分のうちzero-log 5H以外は一般退化因子が被約であり、MPS成分はSNC点を持つため、
この閉条件に全成分は入らない。zero-log23はさらに閉条件D=H⁵に入り、ここには既約三次因子を持つ
両源の点が入らない。他次元の成分とは成分次元から相異なる。
追加half/third25成分もこの次元比較で除かれる。一般rank2成分とも閉rank条件で異なる。
この比較は全空間成分性を証明した後に行い、特殊一点で異なる因子型だけを使った比較ではない。

## 8. 証明書と最初の支持試験の射程

- source producer: `code/classification_nil22_image_stable_cubic.py`, `code/classification_nil22_image_remaining_sources.py`
- 最初の全真正二次支持: `CLASSIFICATION_NIL22_IMAGE_STABLE_QUADRATIC_20261107.json`, `CLASSIFICATION_NIL22_IMAGE_REMAINING_QUADRATIC_20261107.json`
- 二成分の全D/F/真正切片/functional/根基恒等式再構成: `code/classification_nil22_image_components.py`, `CLASSIFICATION_NIL22_IMAGE_COMPONENTS_20261107.json`
- 有限契約: `jobs/nil22_image_stable_cubic.json`, `jobs/nil22_image_remaining_cubics.json`, `jobs/nil22_image_remaining_quadratic.json`, `jobs/nil22_image_components.json`

新二成分以外の九源の最初の二次支持はいずれも正次元。double_A,double_B,quadratic_Aでは二次障害自体が零、
double_C,double_D,triple_C,quadratic_Bでも残normal支持は正次元である。
これは実変形の存在でも既知成分への帰属でもない。高次jetや完全Moraへの自動移行はせず、
次は全源を保つ実幾何変形・既知族認識へ切り替える。

stable_three_linesについて一つの固定T01衝突

    (h,h+εq1,x+3h+εq2,y−3h,x−y−3h)

を全六periodの任意regular Taylor係数0,1,2（18未知数）で検査した。
offDではperiodは括弧の値を二つの線形式で割って得られるためregularであり、逆frameのpoleは1、
wedgeのpoleは2なので、この有限次数で負冪消去と極限値は完全に決まる。
全係数系のrank15に対し、定数targetのaugmented rankは16。
最短障害は極限の `{q1,x}` のx²係数であり、この固定衝突では0、targetでは−1/3。
追加の二次磁気m3,m4,m5にも障害がある。
`CLASSIFICATION_NIL22_IMAGE_TORIC_CONFLUENCE_20261107.json` に全係数左functionalを保存した。
これはこの一つの線形超平面pathの非適合であり、T01又はMPS非所属ではない。

## 9. 二つの全九変数源をhalf/third25へ送る実Poisson曲線

ここでは六磁気m0,…,m5と追加n,c,dを一切特殊化しない。
前節の正次元二次coneを延長せず、全Jacobiを満たす実曲線を構成する。

### quadratic_Bの全源

`K=cx+dy+n`、`R=m0+m1x+m2y+m3x²+m4xy+m5y²` と置き、
座標(q1,q2,x,y)で全非零括弧を

    {q1,q2}=R+q1(cx+2dy+n+y/2)−d q2x+ε(d/2+1/4)q1²,
    {q1,x}=−x²/2,            {q1,y}=xK,
    {q2,x}=−x(y+εq1)/2,      {q2,y}=(y+εq1)K−1

とする。ε=0は§5のquadratic_B全源であり、全εで次数高々二、
全Jacobiが零、Pfaffianは−x²/2、modular場は
`x∂q1+(y+εq1)∂q2` である。

ε≠0では旧nil31-image座標(X,Y,Z,W)から

    (q1,q2,x,y)=(Y−W/ε, εX, Z, W),     det=−ε

という実GL4変換を用いる。元の全15源に次を代入する。

    n6=n7=n9=n11=0, n10=−1/2,
    n0=−(εm0+1)/ε², n1=−m3/ε, n2=(c−εm4)/ε²,
    n3=(6d−4εm5+1)/(4ε²), n4=−m1/ε, n5=−1/ε,
    n8=−(εm2−n)/ε², n12=n, n13=c, n14=d.

変換後の**全六括弧**が上記のregular族と一致する。
n5≠0なので、ε≠0ではaccepted nil31-image half25の全十源に実際に入る。
各固定(m,n,c,d)に対してhalf25の閉性を適用すれば、ε=0も同成分に入る。
従ってquadratic_Bの全九変数源を被覆した。

### triple_Cの全源

まず元の五座標(h,q1,q2,x,y)を(x,q1,q2,h,y)へ入れ替え、
新hをCasimirとするEuler補正を行う。これは唯一の三重超平面old_x=0を無限遠へ取る実PGL変更。
補正ベクトルと全恒等式は証明書に保存した。
この新chartで

    K=c+nx+(d+1/3)y,
    R=m0x²+m1x+m2xy+m3+m4y+m5y²

と置き、全非零括弧を

    {q1,q2}=R+q1(c+nx+(2d+1/3)y)−d q2+εd q1²/2,
    {q1,x}=x/3,              {q1,y}=K,
    {q2,x}=x(y+εq1)/3,       {q2,y}=(y+εq1)K−3x

とする。ε=0で変更後のtriple_C全源に一致する。
全Jacobiは恒等零、Pfaffian=x²、modular場は
`−(2/3)[∂q1+(y+εq1)∂q2]`。

ε≠0ではaccepted third25のnormalized座標(a,b,c_old,d_old)から

    (q1,q2,x,y)=(c_old−d_old/ε, εb, a, d_old),     det=−ε

とし、その全十parameter源へ

    A0=−m3/ε, A1=−(εm1+3)/ε², A2=(c−εm4)/ε²,
    A3=−m0/ε, A4=−(εm2−n)/ε²,
    A5=(9d−6εm5+2)/(6ε²),
    u=n, v=c, t=d, q=−3/ε

を代入する。全六括弧が上記regular族に一致する。
q≠0なのでε≠0はaccepted third25のrank4全源に入り、閉性からε=0も入る。
最後に実PGL/Euler恒等式を戻せば、元のtriple_C全九変数源がthird25へ被覆される。

共通focused checker `code/classification_nil22_image_component_limits.py` は二つの
実GL/source恒等式、全九変数の端点、degree≤2、全strict Jacobiと全350 projective Jacobi、
実Pfaffian、全modular場を直接検算した（0.786秒）。
証明書は `CLASSIFICATION_NIL22_IMAGE_COMPONENT_LIMITS_20261107.json`、契約は
`jobs/nil22_image_nil31_half_limit.json`。この二被覆は成分を追加しない。
別作者reviewは別途行い、統合担当はそのPASS後に受理する。

### double_Dの全源

全六磁気のRをquadratic_Bと同じとし、a,uを§3の追加二変数とする。
次の族は全εでPoissonであり、ε=0でdouble_Dの全源に一致する。

    {q1,q2}=R−[1+(a+1/2)u+(3a+1)y]q1−εa q1²/2,
    {q1,x}=a x²,
    {q1,y}=−(2a+1)x(y+u/2),
    {q2,x}=x+a xy+εa q1x,
    {q2,y}=−u−2y−(2a+1)y(y+u/2)−ε(2a+1)(y+u/2)q1.

Pfaffianは `−x²(y+u/2)`、modular場は `x∂q1+(y+εq1)∂q2`。
ε≠0でquadratic_Bと同じ実GL4変換を用い、元nil31-image15源に

    n0=−m0/ε−u/ε², n1=−m3/ε, n2=−m4/ε,
    n3=−m5/ε−(5a/2+1)/ε², n4=−m1/ε,
    n5=−u/ε, n6=0, n7=−2/ε,
    n8=−m2/ε−[3+(a+1/2)u]/ε²,
    n9=n11=0, n10=a,
    n12=−(2a+1)u/2, n13=0, n14=−1−2a

を代入すると、全六括弧が上記regular族になる。
これは `MODULAR_NIL31_COMPONENT.md` §1の全n7開のcompatibilityを全て満たす。
n7≠0なのでε≠0は既済CN31の25次元成分に入り、その閉性から全八変数のε=0も入る。
a,u又は磁気parameterの非零条件は一切使わない。
`code/classification_nil22_image_double_d_limit.py` と
`CLASSIFICATION_NIL22_IMAGE_DOUBLE_D_LIMIT_20261107.json` は全source/GL恒等式、全Jacobi/Pf、
端点と次数上界を0.458秒で検証した。契約は `jobs/nil22_image_double_d_limit.json`。
この被覆も新成分を追加せず、別作者review後に受理する。

残る未被覆六源はstable二源、double_A/B/C、quadratic_Aである。

## 10. quadratic_Aの磁気正規化と二つの不適合なnormalcone

quadratic_Aについて `α=λ−1/2,β=−λ−1/2` と書く。
全六磁気式F0に対し、fibertranslation
`q1'=q1+A+Bx+Cy`,`q2'=q2+D+Ex+Fy` が加える二次式は、
係数順(1,x,y,x²,xy,y²)で正確に

    (F−B, αD, −βA, 2αE, −(F−B), −2βC).

従ってαβ≠0なら

    A=m2/β, B=0, C=m5/(2β),
    D=−m1/α, E=−m3/(2α), F=m4

によってF0を定数I=m0+m4へ変換する。I≠0ではfiber共通scaleκ=1/Iと
tensor scale1/κによりF0=1となる。これは実affine/tensor同値である。
従って全既約七変数源の稠密開は一変数λの標準源F0=1で尽くされる。
将来この全λ稠密開の既知有限閉和への帰属が得られれば、αβI=0を含む全源へ閉性で延びる。
この正規化だけから帰属は従わない。

既存T36の固定period/smoothing=1 tensorに対し、以前に全Schubert cellで得た
Lagrangian-planeのうち、制限五次式がH³rank3Qとなる二つのprime曲線を試した。
全original座標を(x0,…,x4)、plane座標を(h,x,y)、t≠0とする。

    013prime1:
      (x0,x1,x2,x3,x4)=(x,y,3t^4*x/2000+t*h,h,−3t²*x/100−10h/t),
    014prime2:
      (x0,x1,x2,x3,x4)=(x,10t²h−15t^4x,y,h−3t²x,t*h).

二つのnormalをε倍、tensorをε倍した実normalconeは正則だが、
unique tripleHを∞とするcanonical affine modular線形部は両方rank1、A²=0、
定数部はimage外である。具体的にはfiberを(z,w)とすると

    013prime1: N=(3t^4*x/2500+t/5)∂z+(−3t²*x/250+2/t)∂w,
    014prime2: N=−4t²∂z+(−6t³*x/5−t/5)∂w.

これはquadratic_Aのrank2 nil22-image型と異なる。
残磁気項を作るためのplane-motionはこのnormalanchor/cofactorを変えないので、
この二候補は磁気生成の段階に入る前に不適合と分かる。
同じnormalcone探索を別の点・T番号へ増殖せず、この経路を止める。
`jobs/nil22_image_quadric_normal_cone.json` に元tensor・二plane式・射程を保存した。
これは全MPS非所属の主張ではない。
