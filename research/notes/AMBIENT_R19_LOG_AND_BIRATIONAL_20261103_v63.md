# Rの対数型成分への包含・有理Casimirs・二次CremonaによるPoisson有理同値

2026年11月3日。全Rは既存の四つのPoisson成分CA28,CB28,LPB39,Log(1,1,1,2)28に入る。
以下は新しい成分を加える結果ではなく、既存成分間の交叉の認識である。

## 定理

1. **全R19⊂Log(1,1,1,2)28** である。LPB39への既存包含と合わせ、全Rの各点は少なくとも四つの既知Poisson成分を通る。
2. R\Nの全Poisson構造は、実際の二次Cremona変換によりQ=0の標準形へPoisson有理同値となる。
   一変数のprojective modulus κ²は、この有理同値では消える。
3. R\Nの全標準形は二つの独立な有理Casimirsを持ち、一般葉のZariski閉包は有理曲面である。
   一般のAF≠0標準形では、一般葉閉包はL=C∂yをvertexとする次数6の有理曲面coneである。

有理同値と葉の記述をR∩Nの全特殊点へ広げない。
全Rへの対数型包含は、閉性により別に全境界まで成立する。

## 1. 標準形と一つのDarboux二次式

全R\Nの標準形を

    P_Q=∂y∧V_Q,
    V_Q=4x²∂x−2xz∂z+(2xw+Q(h,z))∂w,
    Q=A h²+Bhz+Fz²

とする。

    K2=120xw+20A h²+15Bhz+12Fz²

と置けば、正確にV_Q(K2)=6x*K2である。
従って二つの同次次数零の有理函数

    phi=x*z²/h³,
    psi=z*K2/(120*x*h²)

はいずれもV_Qで零化され、yにも依存しない。
よってP_Qの有理Casimirsである。

## 2. 旧Log1112 familyとの全係数一致

divisorsを(h,x,z,K2)、degree vectorを(1,1,1,2)とする。
二つの閉logarithmic formsを

    alpha=dlog(phi),  residues=(-3,1,2,0),
    beta=dlog(psi),   residues=(-2,-1,1,1)

と取る。両residue vectorはdegree vectorに直交し、Euler-horizontalである。
そのwedge係数は

    (B01,B02,B03,B12,B13,B23)=(5,1,-3,3,1,2).

これは旧 `legacy/v14/code/rank2_lower_bounds/verify_rank2_components.py` の
log_bivector familyで(a,b,c)=(5/4,1/4,3/4)とした全係数に正確に等しい。

Omega_PをEul∧Pのvolume-dual二形式とすると、全十係数について厳密に

    h*x*z*K2 * alpha∧beta = −60*Omega_P

が成立する。
さらに旧builderが使うd(omega)のdual canonical bivectorは、全150係数で

    P_log=−300*P_Q,canonical

となる。同じfoliationを持つという点集合の議論だけではなく、実際のPoisson tensorの定数倍一致である。

quadricが退化したりdivisorsが重なる場合も、log familyの全多項式表示の特殊化としてその閉包に入る。
Rのt≠0 sourceは全てP_Qへ移り、これはRで稠密である。
Log1112成分は閉でPGL5不変なので、その全飽和閉包Rが含まれる。

旧Log1112が全Poisson空間の実際の28次元成分であることは、既存のtangent/family affine29 certificateを入力として使う。
その大きな計算はやり直さず、sourceとcertificateをhashしている。

## 3. 二次Cremonaは実際にPoisson tensorをQ=0へ送る

次の有理座標を取る。

    w_tilde=w+A h²/(6x)+Bhz/(8x)+Fz²/(10x)=K2/(120x).

h,x,y,zは固定する。この変換の逆は同じ三項を引くだけである。
正確に

    V_Q(w_tilde)=2x*w_tilde

となり、他の成分は変わらない。
従って全五変数のJacobianを使ったactual pushforwardは

    Phi_*P_Q=∂y∧(4x²∂x−2xz∂z+2x*w_tilde∂_(w_tilde))=P_0

に厳密に等しい。
Phiは同次次数1の有理写像としてEulerを保つためP4へ降り、homogeneous quadratic表示は

    [h:x:y:z:w] ↦ [120hx:120x²:120xy:120xz:K2]

である。Q≠0ならdegree2、Q=0ではcommon factorを除いて恒等写像となる。
projective scalarの通常正規化はyのlinear scaleでも吸収できるので、これはPoisson有理同値である。

## 4. 一般葉の有理parameterization

h=1 chartでu=z/h、Y=y/hを取り、二つのCasimirの値をr=phi,t=psiとする。
逆式は

    x/h=r/u²,
    z/h=u,
    w/h=t/u−A u²/(6r)−B u³/(8r)−F u⁴/(10r),
    y/h=Y.

である。r,t,u,Yが全函数体の有理座標となるので、generic leafはu,Yの二変数で有理parameterizeされる。

F≠0の一般葉を商P3へ射影し、上の四座標をu²で払うと

    [u²:r:u³:t*u−A u⁴/(6r)−B u⁵/(8r)−F u⁶/(10r)]

となる。r≠0,F≠0ならP1への同次化にbasepointはなく、degree6である。
u=z/hを戻せるため像へbirationalで、商曲線は次数6の有理曲線である。
元のP4の葉閉包はそのlinear-projection inverse、すなわちvertex Lを持つ次数6の有理coneとなる。

## 証明書と範囲

* `code/ambient_R19_log_and_birational_20261103_v63.py`
* `reports/AMBIENT_R19_LOG_AND_BIRATIONAL_20261103_v63/report.json`

logarithmic十係数、domegaの全canonical150係数、二Casimirs、Cremonaの全Jacobian pushforward、
Euler保存、函数体の全逆式をQQ(A,B,F)上0.23秒で確認した。
有理同値・葉の結論は全R\Nに対するもので、R∩Nの全境界の軌道や全Poisson分類を尽くすものではない。
