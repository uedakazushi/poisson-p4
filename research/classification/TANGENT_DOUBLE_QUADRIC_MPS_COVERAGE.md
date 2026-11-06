# 階数4二重二次曲面の接触配置を覆う二つの実変形

2026-11-07（Asia/Tokyo）。`DOUBLE_QUADRIC_SUPPORT.md` 第7節で得た
isotropic の一変数族の閉包と anisotropic の一軌道を、既知MPS成分へ接続する。
接ベクトルだけで所属を推測せず、全パラメータでPoissonとなる実際の族を使う。
本ノートは一般ファイバーの五点尽くしや、他の退化因子型の分類を主張しない。

## 1. Isotropic 一変数族の全体はT18に入る

座標を(h,u,v,z,w)、hをCasimirとし、c≠1、lambda=2c−1、mu=1−c、rho=2−2cとして

    P_c,uv=lambda*h*u,          P_c,uz=−mu*h*w−c*u*v,
    P_c,uw=−lambda*u²,          P_c,vz=−mu*v²−rho*h*z,
    P_c,vw=−h*w−rho*u*v,        P_c,zw=c*v*w−rho*u*z

とする。既に全cについてJacobiとPf(P_c)=(c−1)(hw+uv)²が確認されている。
以下P_2をc=2の点とする。

### 全次数の実Poisson pencil

hをCasimirとする二次bivector Rを、同じ六成分順で

    R=(-uz−vw/4, −wz/4, 3w²/8, −z²/2, wz, 0)

と定める。直接多項式計算で[P_2,R]=0、[R,R]=0であり、従って全tでP_2+tRはPoissonである。
そのPfaffianは正確に

    Pf(P_2+tR)=(hw+uv−t*wz/4)²−(t/8)(2uz−vw)².

t=8とし、線形座標

    X=u−w,  Y=u+w,  Z=v−2z+h/2,  W=v+2z−h/2

へ移ると、行列の全25成分が次のh-Casimir liftに一致する。

    {X,Y}=−6XY,       {X,Z}=6XZ,        {X,W}=−2XW+2hY,
    {Y,Z}=2YZ+2hX,   {Y,W}=−6YW,      {Z,W}=4ZW+h².

さらに頂点(h,X,Y,Z,W)を、旧番号から新番号への置換

    (0,1,2,3,4) -> (4,1,3,2,0)

で並べ替えると、既存の**完成した**T18 sourceの

    (p0,p1,t0,t1)=(10,12,−2,−2)

におけるbivectorに、Euler wedge diag(0,2,−4,6,0)xを加えたものと完全一致する。
後者は射影bivectorとして零である。入力は
`reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T18.json`。
period Pfaffianと補正分母2p0−p1はいずれも8で非零なので、源の未定義点を使っていない。
この同一視は、平滑化前のトーリック点とT18の一般構造を取り違えたものではない。

対角一径数群の重み(h,u,v,z,w)=(0,1,−1,−2,0)を使うと、pushforwardの規約で

    g_lambda* P_2=lambda^(-1)P_2,   g_lambda* R=lambda R.

従ってlambda≠0でlambda*g_lambda*(P_2+8R)=P_2+8lambda²R、lambda=0の極限はP_2である。
T18は閉じたPGL不変な成分なので、P_2もT18に属する。

### 一点の所属から全cへの正当な接続

P_2の**全150係数**による射影Jacobi微分について、独立実装が350×150行列を再構成し、
102角minorを法1000003で471297≠0と確認した。全係数は有理整数である。
従って標数零の射影接空間次元は125−102=23以下。
既知T18は23次元の全空間成分でP_2を含むため、局所次元も23以上である。
局所次元と埋込み次元が一致し、P_2は全Poissonスキームの正則点となる。
特に、P_2を通る既約成分はT18ただ一つである。

c↦[P_c]の像の閉包は既約でP_2を含む。有限個の全空間成分のうち一つに含まれ、
P_2での唯一性からその成分はT18である。従って**全c≠1**がT18に属する。
cの共鳴値を一般点と同一視した議論ではない。

`DOUBLE_QUADRIC_SUPPORT.md`のisotropic全解の既約五変数源は、一般点でこの族のPGL飽和へ入る。
その全像も閉性によりT18に入るので、isotropic枝の全支持が被覆される。

## 2. Anisotropic 一軌道も既知MPS成分の和に入る

同じ座標でh-Casimir lift Pを

    P_uv=hu,                 P_uz=−hw−3u²−v²,
    P_uw=−2uv,              P_vz=−2hz−2uv,
    P_vw=−3hw−u²−3v²,       P_zw=4uw−4vz

とする。Pf(P)=−3(hw+u²+v²)²であり、これが第7節の唯一のanisotropic軌道の代表である。
ここでは、次のRはhをCasimirとしない。全150係数の射影条件を保持する。
番号0,1,2,3,4はh,u,v,z,wに対応し、非零成分を

    R_01=−hu/2, R_02=−3hv/2, R_03=hz/2, R_04=3v²/2,
    R_12=−uv/2, R_13=uz/4, R_23=vz

と定める。全tについてE wedge[P+tR,P+tR]=0が恒等的に成り立つ。
その退化五次式は

    −3/8 * (tuz+4hw+4u²+4v²)
         * (thvw+tv³+2h²w+2hu²+2hv²).

従ってこれは接空間内だけの候補ではなく、実際の射影Poisson pencilである。

### SNCトーリック点への実PGL極限

P+Rに対して重み(h,u,v,z,w)=(1,0,0,−3,−2)を使う。
二次bivectorの各単項式の重みを、係数単項式の重みから二つの微分添字の重みを引いたものとする。
全単項式の重みは非負であり、重み0の初期bivectorは上のRからR_04だけを除いたものとなる。
その退化五次式は

    −3huvzw/8.

五つの座標超平面は独立でSNCをなし、初期bivectorは一般階数4のlog-canonical構造である。

ここで用いるのは、SNC点への特殊化だけから元の点の所属を推論する誤った論法ではない。
P+Rを含む任意の全空間成分Cを取る。PGL5は連結なのでCを保ち、Cは閉だから
上のトーリック極限もCに属する。MPS Theorem 1.4の「SNC log symplectic locusと交わる
全Poisson成分は既知40個」という主張により、この点を通る全成分は既知40成分の中にある。
従ってCはその一つであり、P+Rは既知MPS成分の和に入る。
この適用の一次資料と仮定は`LOWER_DEGREE_SNC_COVERAGE.md`およびMPS40のcrosswalkと同じである。

次に重み(0,1,1,2,2)を用いる。上の「係数から微分重みを引く」作用rho_lambdaでは

    rho_lambda(P)=lambda^(-1)P,  rho_lambda(R)=R,
    lambda*rho_lambda(P+R)=P+lambda R -> P.

既知MPS成分の有限和はPGL不変かつ閉なので、Pもその和に属する。
第7節で全anisotropic解がこの一軌道へ帰着するため、その全支持が被覆される。

トーリック極限の最大平滑化図式は、置換(2,4,1,3,0)でT33の図式と一致する。
この図式の一致だけを、P+RがT33の完成sourceの一点だという証明には使わない。
ここで必要かつ証明した結論は、既知MPS成分の有限和への包含である。

## 3. 二つの源が全接触配置を覆うことの独立点検

`DOUBLE_QUADRIC_SUPPORT.md`第7節の正規形を、上の族を先に仮定せず確認した。
rankΩ=2の場合、核Kに値を持つaffine場X=Ly+aについてdivX=0、XQ=−kQとなる。
Hessianの核Rは二次元でL不変、商のtraceは−k、R上のtraceはkである。
一次部分のR制限ellが非零でell L=−k ellなので、L|Rの固有値は−k,2k。
k≠0よりR=L(R)⊂imL⊂K、従ってK=Rとなり、商のtrace−kに矛盾する。

Ωがsymplecticの場合、S=Ω^(-1)Hに対する[L,S]=(1−c)SからSはnilpotentとなる。
rankS=2で、symplectic条件はJordan型(3,1)を排除する。
実際、長さ3の鎖a,Sa,S²aの端は鎖の全てと残りの核ベクトルにΩ直交してしまう。
従って型は(2,2)、kerS=imSはLagrangianである。
Ω標準座標ではQ=q2(x)+bᵀy、b≠0とでき、conformal symplectic条件とZQ=cQは
A=cI/2+so(q2)、Ab=(1−c)bを与える。ゆえに(2−3c)q2(b)=0。

q2(b)=0の場合、Q=w+uvに対してZQ=cQとLΩ=Ωの全係数比較を行うと、
第7.3節のc,p,q,r,eによる五変数表示が過不足なく得られる。
一般cでLは可逆かつ単純固有値を持ち、唯一の零点と固有symplectic座標でP_cへ移る。
重みcの二次項に別の項が許される共鳴値は、稠密開から除いてよい。
例えばc=4/5で可能なz²項もHessian rank2から零となる。
既約源の全点はこの一般部分の閉包に入るので、特殊なcも除外されない。

q2(b)≠0の場合、c=2/3、so(q2)の元が非零bを殺すのでA=I/3となる。
下左block Cは対称で、symplectic shear y↦y+Sxにより、固有値差1/3を使って除ける。
Lは可逆なのでaffine定数項を唯一の零点への平行移動で除ける。
二次形式と非等方ベクトルの複素直交標準形により、残る代表は第2節のPの一軌道である。
この独立点検により、上の二つの実変形が単なる選んだ例でなく、全接触配置の支持被覆に使える。

## 4. 独立確認と射程

T18側の独立記録は
`reports/CLASSIFICATION_TANGENT_DOUBLE_QUADRIC_T18_REVIEW_20261107.json`。
全c・全pencilパラメータのJacobi、Pfaffian、全座標変換、完成sourceとの全係数一致、
orbit極限、および全150方向の接空間上界を別実装で確認した。確定版の実測0.458秒。

anisotropic側は
`reports/CLASSIFICATION_ANISOTROPIC_MPS_DEFORMATION_REVIEW_20261107.json`。
疎有理係数の外積と微分により、全パラメータの射影Jacobi、五次式、二つの重みとSNC極限を確認した。
実測0.021秒。どちらもchecker全文・hash・入力hashを保持する。
有限体での零判定を標数零の恒等式へ移しておらず、有限体を使ったのは非零minorの下界だけである。

以上と`DOUBLE_QUADRIC_SUPPORT.md`の全点認識を合わせると、rank4の接触配置も被覆される。
これは指定した二重二次因子層の支持結論であり、新しい成分を追加する結果ではない。
rank3接触配置、他の退化因子型、全S61の一般五点尽くし、完全分類は別の未解決義務である。
