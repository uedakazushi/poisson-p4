# 非対称Q_b側actual D全環の独立承認

2026-11-03。`AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63` をQ1と非空unit open上で承認する。
先行P_a型のlabel交換ではない。新しい三source D3物理監査を明示的に依存させ、全121 ambient odd vectorを保持する。

## 座標とcore81の下包含

CAはsmooth29、CB2との交叉はsmooth28のN=(f)。CB1との交叉Z4は、受理済み相対CI4。
CB二枝meeting C2はdouble smooth Sで、CAに入る：CAの式はSで零、そのC2-nil方向の微分も
全Sで零なので、square-zero ring C2への制限が恒等零となる。
実quotient diagramにより、Z4内のN方程式はisotropic kernel方向である。
従ってCA環Bを `P[[z7,f,rho]]`、P=C[[s20]] と取り、
`J_Z=(z7,f²,rho²)`、C2=Z4/(f)と同時に正規化できる。

CBのminimal38変数presentationはbranch1 `(v8,tau)=0`、branch2 `(u8,tau-rho²)=0`。
Nはbranch2内に入り、CAのN横断方向はbranch2接空間を外れる。
一方T_CAは二CB枝の接空間の和37に入るので、d tauはT_CAを殺す。
従ってCAへ制限したu_iの少なくとも一つはf×unitであり、全u_iとtau-rho²はfの倍数。
全81 core relationの像は正確にf*K、`K=(v8|CA,tau|CA)` となる。

ここでK=J_Zを無条件に仮定しない。
CBの追加common ambient equationsはwhole tangent38を殺すので、CAへの制限にlinear項がない。
従ってJ_Zの七linear pivotはKのもの。
独立コードは、前段の独立source Jacobian・真二次差・独立raw normal quadraticを再利用し、
全150×3がwhole38 blockに入ることをQQ(b)で新たに完全確認する。
これにより追加common equationの二normal quadraticは零。
Z4のfull quadratic rank2はKの二quadricから来て、そのdegree3 spanがunit。
七linearを消去した後に、前段と同じNakayamaでKが全J_Zとなる。
従って全actual D ideal Qには**全idealとして**fJ_Zが入る。
非flat base change後のfiber-product exactnessを使った推論ではない。

## 全ambient上包含と最後のsocle

D_red=N、Z4⊂Dから `fJ_Z subset Q subset (f) intersect J_Z`。
両端の差はP-free rank2、basis f²,f²rho。
全S上のquadratic-zero gateがf²の係数を除くため、残りはfree rank1 P*f²rhoのsubmodule。
Pはdomainなのでgenericに零のsubmoduleは零。

同じ29個の**actual target** baseは三sourceでétale。
独立コードでcached physical Jacobianから三つのprojection detを再計算し、
CB1=`−1728 b8(b²−5)/125`、CB2=`−442368 b9/25`、CA=`8064 b/25`を確認する。
CB ringは同じR上のsplit `R directsum I*t,t²=1` であり、全ambient residual121座標がb_j*tを与える。
CAとの差をf*c_jとすると、c0=1へ可逆row操作後の全ideal recipeは
`f*(J_B,f−2b0,(J_B:b0))`。productsだけでは不十分である。
regular9のKoszul syzygyと、全redundant b_jをその九元へ分解する関係が全syzygyを生成する。
それらのCA像はf倍のc-defectであり、前段の独立三次監査は全150・全121について
目的次数のcofactorを網羅し、全quadratic defectがJ_Zに入ることを確認した。

本REVIEWは、その新しい独立source curve archiveから全150のZ4 quadraticも追加で復元し、
保存されたactual Z4 matrixへ全係数一致させる。元の三次物理計算は繰り返さない。
J_Z quotientはnormal maximal idealの三乗が零なので、それ以上のnormal項は既に消える。
したがって残るf²rho係数は全ambient equationでQQ(b)上恒等零。
intrinsic G-equivarianceとdominant G×b像からPのgeneric点でも零となり、上のtorsion-free性で
残余全体が零。従って `Q=fJ_Z` が全completed idealとして成立する。

## 全環とprimary・nil・depth

actual ringは

    P[[z7,f,rho]]/(f*z7,f³,f*rho²)。

monomial LCMによって `(f*z7,f³,f*rho²)=(f) intersect (z7,f³,rho²)` を独立に確認する。
後者はS-primaryで、そのgeneric Artin長は3×2=6。
両成分はirredundantなのでAssは正確にminimal Nとembedded S_Q。
nilradical(f)は、fが正則なB内でcolonを取って `B/J_Z` に同型。
P-basisはf,f²,f*rho,f²rhoの四つ、nil²はf²,f²rhoの二つ、nil³=0。
nil moduleの長4をprimary quotientの長6と混同しない。
socle f²rhoは非零で全九normal変数に殺され、core depth0。Pの20正則因子からdepth20。
Hilbert seriesは `1/(1-T)^28+(T+2T²+T³)/(1-T)^20`、dim28/edim29/e1。
projectiveはdim27/edim28/depth19/e1。

## 大域embedded閉包とP_a側との区別

S_Qはactual G×b像のirreducibleな代数的閉包。
完成環のembedded primeは、このalgebraic smooth supportの完成idealなので、そのcontractionは
代数的local actual Dのassociated primeとなる（忠実平坦completionのAss公式）。
localizationを戻してglobal associated component S_Q20/projective19を得る。
P_a側はgeneric whole CA/CB旗数(2,1)、Q_b側は(1,2)で、双方に独立の全旗監査がある。
もし二つのirreducible閉包が同じなら、両constructible dense imageは共通の非空openで交わり、
whole旗数が矛盾する。従って二つのembedded prime closureは異なる。
全embedded primeのexhaustionや新しいPoisson成分を主張しない。

範囲はQ1と検証済みunit openおよびそのG-saturationのみ。
b=0/B18原点、全境界点、元のP_a型のlabel交換、global primary listの完了は含めない。
独立コードとhash証拠は `AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW` に保存する。
