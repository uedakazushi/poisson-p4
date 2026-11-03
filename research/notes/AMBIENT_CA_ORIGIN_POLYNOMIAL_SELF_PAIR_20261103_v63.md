# B18原点を通るCAの実多項式self-pair

2026-11-03。四次jetで得たt=−1の有限witnessは、実際には全次数で成立する有限多項式のpaired-sourceへ延長できた。
これは**個別CAの異なる二つの旗**であり、CA/CB間のR19 pairを転用したものではない。
CA原点の全局所環・全conductor・全枝を分類したとは主張しない。

## 実多項式の構成

座標を(h,x,y,z,w)、M=E10−2E03とする。M³=0なので

    g(a)=I+aM+a²M²=(I−aM)^−1,
    g(a)^−1=I−aM

はQQ[a]上の実際の互いに逆な5×5行列である。det g=1。
独立監査済みCA fixed-flag行列C150×20を用いて

    P(a)=P0+a*(-3*C[1]−C[9])

と置く。Cの添字は零始まり。これは四次jetのt=−1,y=0の第一sourceと同じである。
物理テンソルとしては

    P(a)=∂y ∧ V(a),
    V(a)=(x²+12a*h*x/5)∂h−3a*x²∂x/5
           +(h*x+2a*x*z/5)∂z+(h²+z²−8a*w*x/5)∂w.

全150係数のg_*Pを、外積行列gと逆座標g^-1の二次monomial代入から直接計算した。
Lie作用、logarithm、有限jet打切りはこの確認には使わない。
計算上の完全次数上限は7だが、実際のg_*Pはaの二次式に簡約する。
全150係数で

    g_*P(a)=C*second20(a)

がQQ[a]上で恒等的に成立する。second20もcertificateに保存した実多項式である。
従って第一・第二sourceは両方とも同じ実CA fixed-flag fiberに属する。
両テンソルの25 divergence係数と350 strict homogeneous cubic Jacobi係数も、全QQ[a]多項式として零である。

## 同じP(a)に属する二つの異なる旗

critical packetの基準旗は

    F0: L=<e_y>, W=<h,x>⊂ann(L), H=<x>⊂W.

P(a)はF0で表示され、g_*P(a)もF0で表示されるため、同じP(a)はF0とg^-1F0で表示される。
後者の旗は

    L'=<e_y>,
    W'=span(h−2a*z, x+a*h−2a²*z),
    H'=span(x+a*h−2a²*z)

である。a≠0ならH'のh係数がaなのでH'≠<x>、従って二つの実旗は異なる。
これは全てのa≠0に対する明示的な旗であり、単なる接方向の差ではない。

P0とP1の物理150係数のrank2 minorはrows(16,20)で−12/5である。
さらにaに依存しない非零の物理係数も保存した。
従ってP(a)はa=0を通る非定数・非零のordinary tensor arcである。

## 結論と限定

この構成は実際のQQ[a] paired-sourceであり、有限jetの可解性からの推測ではない。
任意次数までのordinary self-pairの不存在という主張は、このarcによって否定される。
ただし、この一本のarcからCA原点の全odd ideal、全conductor、全枝の個数・接合型は決まらない。

この族はconstant vector factor∂yを持つ。従ってgeneric Poisson rank4のself-pairを構成したものではない。
CAのfixed g=C[3]係数は恒等的に零である。
generic source immersion、全proper fiber、cross-tangent、およびこのrank2族の幾何的位置は別packetで検討する。
Nの異なるsource模型にこのCAの旗を無検証で同一視しない。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/ambient_CA_origin_polynomial_self_pair_20261103_v63.py --output <new-output-directory>

sourceと全150物理恒等式:
`code/ambient_CA_origin_polynomial_self_pair_20261103_v63.py`、
`reports/AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/`。
実行は約0.1秒、制限60秒・2.2GB。GBは使っていない。
