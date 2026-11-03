# 実period源の陰関数証明：全Qaを通るMPS成分の追加

2026-11-03。全有限複素aについて、[Qa]は既存MPS成分
**T04–T15、T18–T21、T28**に属する。
既に保存したT01、T02、T03、T16、T17、CA、CBと合わせ、
全Qaを少なくとも24個の相異なる既存成分が通る。
所属の十分な証明のみで、通過成分・局所枝の尽くしや成分下界57の増加はない。

## 1. 実Poisson arcとperiod源

[凍結済みの追加二項arc](Q_MPS_EXTENDED_TORIC_ARCS_20261103_v63.md)を用いる。
今回は対角固有値を d=(0,1,z,x,y) とする。

    s=3(rho−1)/2, a=8rho/[9(rho²−1)], mu=−4/[9s²(s+3)],
    Nt=Ns+t diag(d), Wk=Nt^k−tr(Nt^k)I/5,
    Pi_t=mu*(alpha*t² W1∧W2+beta*t W1∧W3+rho W1∧W4+W2∧W3)。

可換tracefree行列の外積なので全tでcanonical Poissonであり、
追加二項がt=0で消えるためPi_0=Qaである。
この恒等式は固有値に依存しない。相異なる固有値、t≠0の上では同じ明示的固有frameで対角化できる。
共通非零scalar mu*t^5 を除く六つのbivector係数は

    gamma_ij=alpha(di*dj²−di²*dj)+beta(di*dj³−di³*dj)
             +rho(di*dj⁴−di⁴*dj)+(di²*dj³−di³*dj²)。

そのPfaffianはrho倍の四固有値Vandermondeである。
標識置換後の**period**は線形Pfaffian双対とEuler chart変更によりM_sigma*gammaで与える。
bivector係数へbiresidue条件を直接課してはいない。

## 2. 実period条件の滑らかな有理点

対象Tiの保存period basisを列行列B_i、その完全な線形annihilatorをL_iとする。
六変数(alpha,beta,z,x,y,rho)の実代数的sourceにおける条件は

    F_i=L_i*M_sigma*gamma=0。

各対象についてz=2,x=3,y=4の有理点を発見した。
全F_iの厳密な消失、実period座標B_i*c=M_sigma*gamma、
rho(rho²−1)の非消失と固有値の相異性をQQ上で検査した。
さらにrho以外の五変数に関するJacobianの正方小行列式は次の通り。

|対象|rho|alpha|beta|選ぶ非rho変数|小行列式|
|---|---:|---:|---:|---|---:|
|T04|−10/3|−263/3|85/3|alpha,beta,z,x|−201984|
|T05|1/2|49/2|−8|alpha,beta,z,x|−67392|
|T06|−2|−45|15|alpha,beta,z,x|−117504|
|T07|−3/4|−41/4|7/2|alpha,beta,z,x|26784|
|T08|−5/4|−119/4|55/6|alpha,beta,z,x|22752|
|T09|5/6|74/3|−59/6|alpha,beta,z,x|−131136|
|T10|15/4|277/4|−63/2|alpha,beta,z,x|−130032|
|T11|−15/16|−281/16|45/8|alpha,beta,z,x|23922|
|T12|5/3|109/3|−47/3|alpha,beta,z,x|−189696|
|T13|−5/4|−79/4|15/2|alpha,beta,z,x|28512|
|T14|−3/2|−55/2|10|alpha,beta,z,x|28512|
|T15|5/12|199/12|−37/6|alpha,beta,z,x|−35760|
|T18|5/6|95/3|−32/3|alpha,beta,z,x|−94944|
|T19|−1/8|117/8|−13/4|alpha,beta,z,x|27360|
|T20|1/2|49/2|−8|alpha,beta,z,x|−684288|
|T21|5/3|139/3|−52/3|alpha,beta,z,x|147456|
|T28|−3/4|−41/4|7/2|alpha,beta,z,x,y|369360|

標識置換と全方程式はJSON証明書に保存した。
T04–T21は四つの独立period条件、T28は五つである。

## 3. 陰関数定理が実際の曲線を与える

表で選ばなかった非rho変数はその有理値に固定する。
選んだr変数とrhoのA^(r+1)において、F_i=0のr方程式の選択Jacobianは可逆である。
従ってその有理点でrhoへの射はétaleであり、sourceは滑らかな曲線になる。
その点を通る既約曲線成分のrho-imageは稠密である。
これは**実period条件を定義する有限多項式source**への代数的陰関数定理である。
Poisson tangent vectorが障害なく積分するという仮定ではない。
producerはrhoの一次・二次jetもQQで検査するが、全次数の存在根拠は可逆Jacobianである。

この曲線上の非空開で固有値は相異なり、rho(rho²−1)も非零。
全period条件が正確に成り立つので、各Pi_tのt≠0部分は
対象MPS成分の**実toric-boundary subfamily**にPGL5で移る。
このsubfamilyの所属は保存MPS polynomial familyとその閉包を用いる。
edge-smoothingの分母が消える特別periodでも、純toric polynomial subfamilyの閉包所属は保たれる。
対象成分は閉なのでPi_0=Qaも属する。

rhoが曲線上非定数だからa=8rho/[9(rho²−1)]も非定数であり、a-imageはA1で稠密。
Qaは全有限aで非零の多項式写像だから、閉成分の逆像が閉であることにより
**a=0を含む全有限複素a**へ所属が延びる。
有限体の点や数値近似からの外挿は使っていない。

## 4. 保存と射程

producer `code/q_mps_implicit_toric_arcs_20261103_v63.py`、
係数 `code/data/q_mps_implicit_toric_arc_witnesses_20261103_v63.json`、
receipt `reports/Q_MPS_IMPLICIT_TORIC_ARCS_20261103_v63/`。
全17点、実period所属、厳密Jacobian、非零小行列式、一次・二次rho jetをQQで検査する。
実arcのcanonical性・極限・固有frameの共通証明は凍結済み先行packetを用いる。

他のMPS成分への非所属、全通過成分数、正規化fiber、全Poisson成分の尽くしは未結論。
ここでT16,T17等はMPS40の内部IDであり、R19の別の境界層名とは区別する。
