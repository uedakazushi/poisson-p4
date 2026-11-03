# unconstrained positive-weight6族の幾何：3H+Q、指数3 modular flow、二parameter normal form

2026年11月3日。暫定ID `R4_POSITIVE_WEIGHT6_25`。
本稿は既存99-variable finite graphから**N²条件を課さず**に得る三parameter族の幾何を決める。
全Poisson空間のcomponent認識と新しい個数下界は、独立のambient上界・既知成分排除のpacketに依存し、本稿単独では主張しない。

## 1. 全三parameterがPoissonであること

凍結99 graphの最初の90係数を取り、元のpositive affine tensorへ足す。
center/pencilの九係数は全て0のままである。
旧incidenceの唯一の残余4*n0²*n2/3はN²条件に属する行であり、Poisson条件ではない。
今回この行を課していない。
得られたP(n0,n1,n2)の四つの全affine Jacobi多項式を、QQ[n0,n1,n2,x,y,z,w]上で厳密に0と確認した。

非零独立成分は次の通りである。

    Pxy = n1/36−256w/3,
    Pxz = n1*n2*y/48+n1*x/48−64*n2*w*y−64*w*x,
    Pxw = 8/9,
    Pyz = n0*n2*y²/2+n0*x*y−4n0*z/3
                  −n1²*w/1536+n1*w²−n1*x²/64+64*w*y−2/9,
    Pyw = 4x/3,
    Pzw = n2*x*y+x²+2y/3.

したがって単なるtangent directionや有限jetではなく、全parameterの実Poisson族である。

## 2. actual quinticは3Hと滑らかなquadricの和

affine Pfaffianは二次式であり、homogeneous座標(h,x,y,z,w)で

    5184*Q = 2304*n0*n2*y²+4608*n0*x*y−6144*n0*h*z
             −3*n1²*h*w+4608*n1*w²−72*n1*x²+96*n1*h*y−1024*h²

とすると、元のquadratic bivectorをhomogenizeした**actual quintic**は、固定非零scalarを除き

    h³*Q

である。Pfaffianのquartic/cubic部分の消失を全係数で確認した。
canonical Euler quotientを取ってもactual quinticは変わらない。

Qのsymmetric matrixのdeterminantは

    1024*n0³*n1*(32*n0+n1*n2)/531441.

従ってn0*n1*(32*n0+n1*n2)≠0上でQは滑らかなirreducible quadricで、H:h=0は唯一のmultiplicity3 hyperplaneである。
Hへの制限の(x,y,w) minorは

    −4*n0*n1*(32*n0+n1*n2)/729.

zがradical方向なので、H上のquadricはrank3である。
点v=[0:0:0:1:0]でgradient Q=(-32*n0/27,0,0,0,0)となる。
よってHはQのvにおける**tangent hyperplane**で、H∩Qはvertex vのquadric coneである。

## 3. modular fieldは指数3のunipotent symmetry

affine divergenceは

    Vmod=(4n0/3)*∂y+(n0*x+n0*n2*y+n1/48)*∂z.

そのhomogeneous5×5 matrix Nは

    N_yh=4n0/3,
    N_zh=n1/48, N_zx=n0, N_zy=n0*n2,
    その他=0.

N³=0で、N²の唯一の非零候補はN²_zh=4*n0²*n2/3である。
n0*n2≠0上でrank列は(2,1,0)、Jordan typeは(3,1,1)である。
従ってgeneric memberに旧square-zero incidenceのN²=0を戻してはならない。

有限flow exp(tN)=I+tN+t²N²/2が、元のhomogeneous raw bivectorを全parameter・全tで厳密に保存する。
さらにN^t*Qmatrix+Qmatrix*N=0で、同じflowはquadric Qのorthogonal automorphismである。
符号規約を明示すると

    P*gradient(Pfaffian) = −Pfaffian*Vmod

も全係数で成り立つ。symplectic openではmodular flowは−log(Pfaffian)のHamiltonian flowである。

n0≠0,n2=0ではN²=0でJordan type(2,2,1)、n0=0,n1≠0ではrank1/index2となる。
これら境界での他成分帰属は本稿で分類しない。

## 4. constant pencilと無限遠center

C=∂y∧∂zは全族とcompatibleである。
Pの全generic constant compatible bivectorsはCのspanだけである。
これは全六係数のcompatibility matrixが中心sample(1,2,3)でrank5、Cが全parameterでkernelに入ることから従う。

さらに

    z -> z−3*lambda/(4*n0)

というaffine translationによってPはP+lambda*Cへ移る。
このPoisson pencilは、n0≠0上ではz-translationsで実現される。

compact centerの全linear matrixもsampleでrank4であり、全parameterで

    [u:vx:vy:vz:vw]=[0:0:0:1:0]

をkernelに持つ。従ってgeneric memberのcenterは一つで、**無限遠**にある。
有限centerを仮定するnormal formへ置き換えてはいない。

## 5. 二parameter normal formとgeneric orbit/stabilizer

n0≠0、a≠0、epsilon=±1について、有限変換

    x'=a*x, y'=a²*y, w'=epsilon*a*w,
    z'=a³*z+(a³−epsilon)/(6*n0),
    P'=epsilon*a^(−2)*pushforward(P)

は厳密に

    (n0,n1,n2) -> (epsilon*n0,epsilon*a*n1,n2/a)

を与える。全六成分でrational identityを確認した。
n1≠0ならepsilon=1,a=1/n1でn1=1にできるので、二parameter

    u=n0, v=n1*n2

のnormal formを得る。
さらに(u,v)と(−u,−v)は上のepsilon=−1変換で等価である。
これを全てのfinite PGL同定の分類とは呼ばない。

sample(1,2,3)でcanonical150係数の全GL25 orbit matrixはrank23、scalarを加えるとrank24である。
一方、非零modular Gaが全generic memberのprojective stabilizerに入るため、projective orbitは常に23以下である。
従ってgeneric orbit次元は正確に23、connected projective stabilizerはmodular flowのGaそのものである。
追加のfinite stabilizerは未分類である。

二parameter normal formとorbit上界23から、この族のPGL saturation次元は25以下。
同じsampleのGL+scalar+三parameter微分matrixにrank26のexact minorがあるため、射影saturation次元は**正確に25**である。
そのためorbit quotientは二次元で、上の二parameter normal formはgeneric finite coverを与える。
u,vが全finite PGL同定を既に除いた完全不変量であるとは主張しない。

## 証明書とscope

* `code/r4_positive_weight6_geometry_20261103_v63.py`
* `reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json`

新規QQ計算1.73秒。
全affineP、Pfaffian/Q matrix、Ga finite identity、parameter変換、pencil/center matrices、
全150 orbit/family matrixとexact rank factorizationを保存した。

この25次元irreducible saturationが全Poisson空間の新componentであるという結論には、
別packetのactual ambient上界25と既知成分からの排除を合わせる。
既知下界の更新と成分ID確定はその独立review後に行う。
