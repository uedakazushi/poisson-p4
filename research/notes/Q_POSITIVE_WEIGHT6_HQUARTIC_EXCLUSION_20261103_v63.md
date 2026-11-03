# 正weight6軸は最後の追加四次成分の全閉包にも入らない

2026年11月3日。対象は既知成分 `R4_H_QUARTIC_ESSENTIAL5` の**全閉包**である。
一般退化因子の本質変数数だけでは5H境界を排除できないため、properな追加データを記録する。
以下の有限全次数incidence modelにより、正軸の非零点はこの成分の境界にも入らない。
全Poisson空間の完成環や新成分・既知成分数の増加を主張しない。

## 1. properな三つの追加データ

既知代表は旧v6の `P_*=h²J+Q_*` であり、h=1で原点がcenter、Jは非退化な
constant compatible bivector、affine modular vector field N=div(P_*)は線形でN²=0、rank2である。
この三条件は全affine座標変更と非零P-scalarで保たれる。
コードは元の全代表式からcenter、compatibility、N²=0、rank2を再検査する。

まずPoisson hyperplane Hをmarkする。ある固定chartでH=(h=0)とすると、
Hを保つprojective quadratic Poisson tensorは、そのaffine開上で一意な次数≤2のbivector Pとなる。
これはraw representativeをP^{0i}=0へ取ることと同じである。
H-tangencyによりh=0でP^{0i}=x_i*ellとなるので、Euler vertical項を二度引けば
P^{0i}=0にできる。逆も明らかで、係数環上の定数線形Koszul関係による同定である。
従って全affine係数は6組×15 monomial=90個で尽くされる。

P=C+L+Qをaffine次数で分ける。center p=v/uが有限であるとは、
平行移動後にlinear jetが零、すなわち `u*L+partial_v Q=0` ということである。
この式は[u:v]について斉次線形なので、centerをP4へcompactifyした**閉incidence**を与える。
u≠0が実際のfinite-center locusである。

同時に非零constant bivector CcをP5にmarkし、[Cc,P]=0を課す。
さらにaffine vector field div(P)を5×5行列N(P)としてN(P)²=0を課す。
これらは係数の多項式閉条件である。
従ってmark H、compact center、projective pencilからなるsourceは、
H-chart上でP4×P5をfiberとするproperなclosed incidenceである。
既知四次成分の稠密orbitはこのsourceのu≠0かつCc非退化の部分から来る。

## 2. 正軸点の全preimageは一つ

eta=1の正軸点P+を使う。全非零etaは既存のweighted PGL orbitにより同一視できる。
actual quinticは非零定数倍のh5なので、そのlinear factorとしてのmarked Hの
**幾何的**候補はH=(h=0)だけである。marked-H fiberの被約性は仮定しない。

全compact center条件の24×5係数行列を原P+から構成するとrank4で、kernelは

    [u:vx:vy:vz:vw]=[0:0:0:1:0]

だけである。center fiberはこの線形条件では一点である。
全constant compatibility条件の20×6係数行列はrank5で、そのkernelは

    C0=partial_y wedge partial_z

の一つのprojective lineだけである。両rankのunit minorを保存する。
従って未検査のcenter chartやpencil chartへ逃げるpreimageはない。
使うchartはvz=1、Cc^{yz}=1であり、残るcenter/pencil入力は計9変数である。

## 3. 99変数・209個の実際の二次方程式

P=P++Rとし、Rの全affine90係数を独立に取る。
centerを[u:c1:c2:1:c4]、pencilをC0+DCとし、DCのyz成分を0に固定する。
入力は90+4+5=99個。原点はP+と上の唯一のcenter/pencilである。

全方程式を次の順で構成する。

1. Jac(P)=0の全4組×35 monomial、140式。
2. uL+partial_v Q=0の全6組×4、24式。
3. N(P)²=0の全5×5、25式。
4. [Cc,P]=0の全4組×5、20式。

合計209式で、P+でN=0なので全て次数≤2、定数項0である。
証明書は全式の全有理係数を、変数index列と分子分母の形で保存する。
元のPoisson式の一部だけを使った代用ではない。

Jacobianは209×99、rank77。GL4の16方向、translations4、P-scalar1の
**実際のsource作用**から19個の独立なorbit方向を取る。
centerとpencilの変換微分も全て含む。残るkernel方向は3個である。
その3列、19orbit列、77pivot座標列を合わせた99×99 frameのdeterminantが非零である。
全77×77 Jacobian minorと99×99 frame determinantをQQで計算・保存する。

affine群作用ではcenterのuはu/denへ移り、denは原点で1となるunitである。
またpencilのyz成分による再正規化も原点でunit。
Poisson、center、compatibility、N²の**全scheme ideal**を保つ有限作用なので、
上のunit frameは完成局所sourceをorbit19と線形transverse sliceの積に分離する。
非被約構造を無視したset-theoreticなorbitの議論ではない。

## 4. 全次数で終わる77個のimplicit graph

三つのnormal変数をn0,n1,n2とする。
77pivot方程式のJacobian minorは非零定数なので、形式的implicit theoremが適用できる。
その解を有限に求めると、linear kernel入力に加える非線形補正は正確に次の三つだけである。

    r17 += n1*n2/48,
    r49 += −n1²/1536,
    r54 += n0*n2/2.

他の全pivot補正は0。ここでr0,...,r89の順は、六つの上三角pairについて
`1,x,y,z,w,x²,xy,xz,xw,y²,yz,yw,z²,zw,w²` の15係数である。
linear kernel map、全99入力mapを証明書へ保存する。

この二次graphを**原209多項式へ全次数で直接代入**すると、
77pivot式は全て厳密に0、全209のうち唯一の非零pullbackはzero-based row179の

    (4/3)*n0²*n2

だけとなる。原式が二次、graphが二次なので代入の最大次数は4であり、
実行は次数打ち切りではなく全多項式係数の一致である。
従ってimplicit theoremの一意性から、原sourceの完成環は厳密に

    C[[orbit19,n0,n1,n2]]/(n0²*n2)

となる。これはsource incidenceの完成環であり、全Poisson局所環と混同しない。

**全9個のcenter/pencil chart座標は、この有限graph上で恒等的に0である。**
特にfinite-center座標uは、全source完成環で厳密に0となる。
quad initialのu²からuの冪零性を推測したのではない。
先行のcenterだけのthrough5計算ではこの結論が得られなかったが、
必要なconstant-pencilとN²を加えた実際の全sourceについて有限graphが成立した。

## 5. 全閉包からの排除

Noetherian局所環から完成への写像はfaithfully flatなので、完成環でu=0なら
元のsource局所環でもu=0である。従ってこのsource点のあるZariski近傍では
finite-center開D(u)は空であり、この点はD(u)の閉包に入らない。

Hの四方向を戻しても結論は変わらない。
marked-H chartを局所PGL4パラメータで固定Hへ移す有限作用があり、
centerがH上にあること（u=0）はその作用で不変である。
従って19個の固定H affine orbit方向とHの4方向を混同せず、全marked sourceでも
uは局所的に0となる。scalarを除くと正軸のprojective orbit dimension22とも整合する。

仮にP+が既知四次成分の全閉包に入れば、properなH/center/pencil sourceにより、
稠密finite-center orbitからのspecializationはP+上の何らかのpreimageへ持ち上がる。
上で幾何的preimageは全て一つと確認済みであり、その点はD(u)の閉包に入らない。
矛盾する。したがってP+はentire `R4_H_QUARTIC_ESSENTIAL5` の外にある。
非零etaの全点も同じPGL orbitなので、全非零正軸を排除できる。

これは一般点の退化因子F4が本質五変数という事実や、finite centerが見つからないという事実だけの議論ではない。
properなcompact source、その全preimage、全209式の有限implicit graphが必要な部分である。

## 保存とscope

* `code/q_positive_weight6_hquartic_exclusion_20261103_v63.py`
* `reports/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63/{report.json,certificate.json.gz}`

全source式、Jacobian、unit minor、orbit/transverse frame、有限graph、209pullback、
全projective fiber検査、元Hquartic代表のsource条件を保存する。
全Poisson空間の枝・新成分・全既知成分の排除は結論しない。
MPS T04,...,T38の未排除成分は別課題として残る。
