# 消える低冪項を加えた実arc：全Q_aはMPS T16,T17にも入る

2026-11-03。**全有限複素aについて、[Q_a]は既存MPS成分T02,T03,T16,T17に属する。**
ここでT16,T17はMPS40カタログの内部component IDで、R19に関する別の境界層の名称ではない。
四成分はいずれも既存の射影24次元成分である。
既知T01、CA、CBと合わせ、全Qaを少なくとも七つの相異なる既存成分が通る。
全局所枝・通過成分の尽くし、成分総数の新しい加算は行わない。

[二固有値曲線によるT02,T03所属](Q_MPS_TORIC_ARC_MEMBERSHIPS_20261103_v63.md)は凍結保存する。
今回は実arcの、極限で消える二項を利用した簡単な有理パラメータ証明と、二つの追加成分への所属を得た。

## 1. Qa極限を保つ追加項

rho=lambda/mu、s=3(rho−1)/2、a=8rho/[9(rho²−1)]とし、N_sは既存のregular nilpotent行列とする。
N_t=N_s+t*diag(0,1,2,3,d)、W_k=N_t^k−tr(N_t^k)I/5を用い、

    Pi_t=mu*(alpha*t² W1∧W2 + beta*t W1∧W3
             +rho W1∧W4 + W2∧W3)

とする。mu=−4/[9s²(s+3)]である。
追加した二項はt=0で消え、全150係数でPi_0=Qaを保つ。
W_kは可換・tracefreeなので全tでcanonical Poisson。
producerはさらに四つの実arcごとにdivergence25・raw Jacobi350を厳密に検査した。

t≠0では既存のunitriangular eigenframeで対角化できる。
W1∧W2の対角係数はt³、W1∧W3はt⁴なので、追加したt²,tによって全項の対数型係数はmu*t⁵となる。
共通scalarを除いた六係数は

    gamma_ij=alpha(di*dj²−di²*dj)+beta(di*dj³−di³*dj)
             +rho(di*dj⁴−di⁴*dj)+(di²*dj³−di³*dj²)。

power-basis上のskew係数は(alpha,beta,rho,1,0,0)でPfaffianはrho。
従ってgammaのPfaffianはrho倍の固有値Vandermondeで、alpha,betaに依らない。
periodは前ノートと同じ線形Pfaffian双対D(gamma)であり、bivector係数へbiresidue条件を直接課していない。

## 2. 四つの有理係数による実period所属

全て恒等label置換でよい。
T02,T03ではd=4を固定し、rhoを自由パラメータとする：

|成分|alpha|beta|
|---|---|---|
|T02|(973rho+635)/35|−(322rho+155)/35|
|T03|(139rho+75)/5|−2(23rho+10)/5|

T16,T17ではdを自由パラメータとする。

    rho16=(4d−11)(4d−1)/[3(2d−5)(2d−1)],
    alpha16=(12d⁴+54d³−139d²−303d+136)/[3(2d−5)(2d−1)],
    beta16=−(28d³+6d²−250d+81)/[3(2d−5)(2d−1)]。

    rho17=−(4d−11)(4d−1)/[12(2d−1)],
    alpha17=−(60d³−139d²−297d+133)/[12(2d−1)],
    beta17=(8d³+6d²−122d+39)/[6(2d−1)]。

保存MPS period basisの全線形条件へ代入した結果は、各一変数有理関数体上で恒等的に零。
実periodパラメータも基底から解き、全六成分が一致することを検査した。
alpha,beta,rhoを解いた三次小行列式はそれぞれ

    T16: −12d(d−3)(d−2)(d−1)(2d−5)(2d−1),
    T17: 48d(d−3)(d−2)(d−1)(2d−1)

であり、非空な開の上で通常の線形方程式を解いている。
T02,T03でalpha,betaを解く二次小行列式はd=4で420,360である。

## 3. 十分性と全aへの閉包

一般パラメータでd_iが相異なり、rho(rho²−1)も非零なので、実arcのt≠0部分は
対象成分の実toric-boundary subfamilyへPGL5で移る。
そのsubfamilyが既存成分に入ることと、対象成分が閉であることからPi_0=Qaも対象成分に入る。
これはperiod条件の消失だけを十分条件としたものではなく、実Poisson arcを同時に構成した議論である。

四場合ともrhoとaのパラメータ微分は非零で、a-imageはA1で稠密。
Qaはa=0でも非零で全A1に定義されるので、閉成分のQa逆像も閉である。
従って除いた有限個の極・衝突値・coverの特別値を含む**全有限複素a**へ所属が延びる。
有限体の標本や数値的近似から複素数へ外挿していない。

## 4. 保存と射程

全150係数のlimitとEuler-vertical恒等式、各arcのraw Jacobi350・divergence25、
unitriangular eigenframeの漸化式、実period係数、非定数性をQQ上で検査した。
各arcの150係数をtの多項式として保存している。

producer: `code/q_mps_extended_toric_arcs_20261103_v63.py`。
係数: `code/data/q_mps_extended_toric_arc_witnesses_20261103_v63.json`。
receipt: `reports/Q_MPS_EXTENDED_TORIC_ARCS_20261103_v63/`。

他のMPS成分への非所属、全Qa局所枝の数、正規化fiber、全Poisson成分の尽くしは結論しない。
既存の成分下界57は変わらない。
