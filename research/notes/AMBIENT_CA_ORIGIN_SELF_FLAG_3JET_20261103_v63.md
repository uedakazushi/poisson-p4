# B18原点のCA self-flag三次jet：必要条件は非空なconic

2026-11-03。Nの19係数模型での不整合をCAへ転用せず、**CAの全20係数とparabolic断面free10**で再計算した。
全9 relative-flag方向を保持すると、ordinary flag separationの三次必要条件は非空である。
これは有限3jetの結果であり、全次数の形式的・代数的pair branchを構成したものではない。

## 実source、gauge、membership

同じgeneric B18点P0と、独立監査済みCA fixed-fiber行列C（150×20）を使う。
parabolic orbit10方向を除くsource free indicesは

    (1,2,3,9,10,14,15,17,18,19)

である（零始まり）。第一sourceは `P(a)=P0+a*P1+a²*P2+...`、P1,P2はこのfree10に属する。
ordinary relative flagの一次kernel方向を

    M=E_(1,0)-2E_(0,3)

とし、E_(1,0)座標をaそのものに正規化する。
残る8位置は(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(3,2),(4,2)であり、

    g(a)=I+aM+a²B+a³C3+...,  B,C3はこの8方向全てに自由

と置く。これらは全9次元opposite-unipotent flag chartである。
N=g-IはN⁴=0を満たすのでdet g=1で、実際に可逆な群座標である。

同じ固定旗に属する第二sourceが存在する条件は、canonical物理150係数の`g_*P(a)`がCのimageに入ることである。
保存した20次単元blockから作るprojection `pi=I-C*C_left*row_selector` を全150行に適用する。
Euler quotientや少数のraw出力だけでmembershipを代用しない。

## Log(I+N)による有界な係数計算

巨大な多変数physical pullbackを展開せず、既に監査された150係数Lie作用Lを使う。
group logarithmは

    log g = aM + a²(B-M²/2)
       + a³(C3-(MB+BM)/2+M³/3) + O(a⁴).

この式とexp(log g)=gの次数2,3行列恒等式を記号のまま確認した。
特に `L_M²`（作用の二乗）と `L_(M²)`（行列M²の作用）を混同しない。
二次membershipは

    pi[L_M P1 + L_B P0 + (L_M²-L_(M²))P0/2]=0

であり、10個のP1係数rhoと8個のB係数に関する150×18のQQ affine linear systemとなる。
rankは15、consistentで、自由度は3である。

自由変数をv0=B_(0,4), v1=B_(1,3), v2=B_(4,2)と取ると

    rho=(v0+3v1/2, 0, -v2/2, v0+5v1/2+4,
         0,0,0,-v0/2,0,0),
    B=v0*E_(0,4)+v1*E_(1,3)+v2*E_(4,2).

三次には高次source P2/free10とC3/flag8が一次で現れ、そのimageは同じ150×18行列である。
全18高次変数をこのimageで消去し、全150行のcokerを取った。独立な三式は

    q1 = 2*v2*(2*v0+3*v1+5)/5,
    q2 = -(4*v0*v1+4*v0+7*v1²+18*v1+8)/2,
    q3 = -v2*(v0+4*v1+5)/5.

## 正確な必要idealと有限witness

F=4v0v1+4v0+7v1²+18v1+8、H=4v0+7v1+11とすると

    v2 = -(H/6)q1 + (2v2/3)q2 -(2H/3)q3,
    F = -2q2.

逆包含もq1,q3がv2の倍であることから従う。従って全三次coker idealは正確に

    (v2,F)

である。tiny QQ三変数GBも同じ二生成を与えた。
v1+1はこのquotientで単元で、t=v1+1とすれば

    v1=t-1,  v0=-7t/4-1+3/(4t),  v2=0,  t!=0.

従って必要3jetのparameter曲線はQQ上の滑らかなGmである。
v=(-2,0,0)を選び、高次18変数も保存QQ inverseで解いて、全150 membershipの次数2,3が零になる有限witnessを保存した。
その第一sourceの追加CA係数C[2]はa²で-1となる。
このwitnessはCA20とN19の違いを実際に使うので、Nの不整合をCAへ持ち込めない。
全witnessでこの係数が非零だとは主張しない。

同じ追加CA係数の高次値は、全高次自由度に依らず

    (v1+2)*(3v0+7v1+5)/2

と決まる。これはここでの有限source条件であり、conductorの式とは同定しない。

## 結論の限定

このpacketは、relative flag separationがa一次で、第一sourceがO(a)で動くordinary chartの必要条件を扱う。
高いflag separation、別のramification order、全次数liftingを尽くさない。
可解な3jetから実branchを推論せず、仮にこのchartが不整合でも全てのhigher-order pairの不存在とは推論しない。
個別CA conductor、ν方向のnil、R19 conductor、B18原点の全odd idealは未決である。

全実行は60秒・2.2GB制限内（約1秒）で完了した。新しい大規模GBや高次jetはない。
code: `code/ambient_CA_origin_self_flag_3jet_20261103_v63.py`。
report、全150係数行列、affine解、cokerと有限witness:
`reports/AMBIENT_CA_ORIGIN_SELF_FLAG_3JET_20261103_v63/`。
