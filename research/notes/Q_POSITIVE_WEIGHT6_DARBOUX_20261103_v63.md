# 正weight6軸全体の多項式Darboux座標と(2,3) Poisson Cremona変換

2026年11月3日。a=8/9の正軸について、全etaで有効な有限式を与える。
PGL同値や既知成分への帰属を意味するものではない。

## 1. affineでの全多項式同型

h=x0=1上の旧物理座標を(x,y,z,w)=(x1,x2,x3,x4)とする。
原Poisson曲線の六つの非零となり得る括弧は

    {x,y}=−(256/3)eta*w, {x,z}=−64eta*x*w, {x,w}=8/9,
    {y,z}=−2/9+64eta*y*w, {y,w}=(4/3)x,
    {z,w}=(2/3)y+x².

次の有限座標変更を取る。

    X=x,
    Y=y−3x²/4+48eta*w²,
    Z=z−3xy/4,
    W=w.

QQ[eta,x,y,z,w]上でJacobian determinantは恒等的に1であり、全六括弧は

    {X,W}=8/9, {Y,Z}=−2/9, 他0

となる。inverseも多項式で、正確に

    x=X,
    y=Y+3X²/4−48eta*W²,
    z=Z+3XY/4+9X³/16−36eta*XW²,
    w=W

である。両方向の全4座標のcompositionとinverse側のPoisson tensor恒等式も検査する。
これは有限jetでも点での接空間の同定でもない。

唯一のconstant translation方向∂zのHamiltonianとして

    H=(9/2)y−(27/8)x²+216eta*w²=(9/2)Y

が取れ、原tensorで `P*dH=∂z` と直接確認できる。

## 2. homogeneous二次mapと三次inverse

全P4上のrational mapは

    F([h:x:y:z:w])
      =[h² : hx : hy−3x²/4+48eta*w² : hz−3xy/4 : hw]

である。標的同次座標を[H:X:Y:Z:W]とするとinverseは

    G([H:X:Y:Z:W])
      =[H³ : XH² : YH²+3X²H/4−48eta*W²H :
        ZH²+3XYH/4+9X³/16−36eta*XW² : WH²].

全QQ[eta]で両compositionは厳密に

    G(F(h,x,y,z,w))=h^5*(h,x,y,z,w),
    F(G(H,X,Y,Z,W))=H^5*(H,X,Y,Z,W)

となる。Fの五座標は次数2、Gの五座標は次数3で、両gcdは1である。
このgcdの結論は全eta特殊化で保たれる。
Fの共通因子はhの因子でなければならないが、Fの第三座標のh=0でのx²係数は常に−3/4。
Gの共通因子はHの因子でなければならないが、第四座標のH=0でのX³係数は常に9/16。
従って全etaで次数の取消しはなく、正確なbidegreeは(2,3)である。

原homogeneous Pの全150係数を読み、targetを

    P_const=H²*((8/9)∂X∧∂W−(2/9)∂Y∧∂Z)

とする。両方向で全5×5 tensor恒等式

    dF*P*dF^T=P_const(F),
    dG*P_const*dG^T=P(G)

を正確に確認する。従ってF,Gはprojective rational Poisson同型を与える。
Euler項を法とした同値だけで済ませず、このraw homogeneous tensorの恒等式まで成立する。

## 3. scope

全etaに対するpolynomial Poisson同型と、全P4の(2,3) Cremona Poisson同型を証明した。
eta=0も式は有効である。非零etaが一つのPGL軌道を作る事実は別packetによる。
Fの次数が2なので、これをPGL5座標変更として扱わない。
既知成分への所属、全Poisson局所環の枝や非平坦性、新成分や既知下界の増加は推論しない。

* `code/q_positive_weight6_darboux_20261103_v63.py`
* `reports/Q_POSITIVE_WEIGHT6_DARBOUX_20261103_v63/{report.json,certificate.json}`

全入力と出力証明書のSHA256をreportへ保存する。
