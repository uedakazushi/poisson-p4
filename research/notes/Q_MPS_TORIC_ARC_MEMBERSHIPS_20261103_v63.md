# 全Q_aは既存MPS成分T02とT03にも入る

2026-11-03。**全ての有限複素aについて、実際の射影Poisson構造[Q_a]はT02およびT03に属する。**
T01・CA・CBへの既知の所属と合わせ、これら五つの既存成分が全Q_aを通る。
新しい全域成分を数えた結果でも、局所枝や通過成分の尽くしでもない。
T02,T03は保存カタログの内部IDで、どちらも既知の射影24次元成分である。

対象は(Q_a)14=a x0²、(Q_a)23=−2x0²/9、(Q_a)24=4x0x1/3、
(Q_a)34=2x0x2/3+x1²のcanonical代表。
a=0でも非零のrank2 tensorなので、A1_aから射影Poisson空間への写像は全域で定義される。

## 1. 同じQa極限を持つ、固有値の自由な実toric arc

[T01所属の証明](Q_T01_MEMBERSHIP_20261103_v63.md)のregular nilpotent N_sを使う。
そのsubdiagonal linkは1,s,(s+3)/2,−9/(s+3)である。

    rho=lambda/mu, s=3(rho−1)/2,
    a=8rho/[9(rho²−1)],
    lambda=−2a/(9s), mu=−4/[9s²(s+3)].

s(s+3)≠0の開でN_sはregular nilpotentで、
lambda(N_s x)∧(N_s^4 x)+mu(N_s²x)∧(N_s³x)=Q_a。

ここでD=diag(0,1,2,x,y)、N_t=N_s+tDとする。
W_k=N_t^k−tr(N_t^k)I/5を用いた

    Pi_t=lambda(W1 xvec)∧(W4 xvec)+mu(W2 xvec)∧(W3 xvec)

はtについて多項式で、t=0で同じQ_aへ至る。
W_kは可換・tracefreeなのでPi_tはcanonical Poissonである。
新producerでは全150係数でQa極限とEuler-wedge差を検査した。
一般固有値の場合のunitriangular eigenframeも、

    Sij=t^(j−i)*(-1)^(i−j)*prod_(r=j)^(i−1)n_r /
        prod_(k=j+1)^i(d_k−d_j),  Sii=1

の各entryの漸化式で厳密に確認した。
t≠0で固有値が相異なると、Pi_tはPGL5変換を除き標準対数型になる。
その六つの**bivector係数**はmu*t^5倍の

    gamma_ij=rho(di*dj^4−di^4*dj)+(di²*dj³−di³*dj²), 1≤i<j≤4。

## 2. bivectorとbiresidueの規約

保存MPSカタログのperiod_basisはb12,b13,b14,b23,b24,b34という**biresidue**座標である。
gammaへその基底条件を直接課していない。
4×4skew行列の線形Pfaffian双対

    D(gamma)=(gamma34,−gamma24,gamma23,gamma14,−gamma13,gamma12)

を使う。rank4上ではD(A)=−Pf(A)A^(−1)で、D²(A)=A。
共通scalarを許す斉次period空間だから、この線形双対で正確に判定できる。
5×5biresidueはB0i=sum_(j=1)^4 Bijで拡張し、行和を零にする。
全120label置換について、chartを取り直したbivectorの双対が
置換したBと置換の符号だけ異なることを、六つの独立basis上で確認した。

target成分T_iのtoric-boundary subfamilyは、Bがそのperiod線形空間に入り、edge係数が零の族である。
一般BはMPSの実変形族のtoric基点であり、線形Pfaffian双対によるbivector族は全Bで多項式。
従って、periodの別の開条件が落ちる点も閉包により同じ既存成分に入る。
これは「必要条件が消えたので成分所属」と推論する手順ではない。
下記の実arcをこの実toric-boundary subfamilyへ持ち上げていることが十分性を与える。

## 3. T02への稠密なQa所属

次のQQ既約曲線f2(x,y)=0を取る：

    f2=(14y³−9y²−3y+2)x²
       +(−6y⁵+2y⁴−2y³−y²+8y−3)x
       +2y⁵−2y⁴−4y³+6y²−3y+1。

この上で

    rho=−(7y²−y−2)/(3y³+6y²+12y−7)

とする。label置換は恒等でよい。
D(gamma)が保存T02 period basisのspanに入ることを、QQ[x,y]/(f2)での全六成分と
実際の四つのperiod係数により確認した。

rhoは非定数であり、aは具体的に

    8(7y²−y−2)(3y³+6y²+12y−7) /
    [9(3y³−y²+13y−5)(3y³+13y²+11y−9)]

となる。非空なdistinct-eigenvalue開を取り、上の実arcを使える。

## 4. T03への稠密なQa所属

u=x+y、v=xyと書き、次の対称多項式をf3とする：

    f3=−u⁵+3u⁴v+u⁴+4u³v+2u³−12u²v²−3u²v−4u²
       −9uv²−4uv+3u+6v³+9v²−1。

QQ既約曲線f3(x,y)=0上で

    rho=−(u²−u−2v²−v)/(u³−2u²v−2uv+2v²−1)

とする。これも恒等label置換でD(gamma)がT03のperiod空間に入る。
全constraintsの剰余と、実period係数の六成分一致をQQで確認した。

## 5. 非定数性・非空な開・閉包の論理

両fはxについてprimitiveである。
T02はy=2、法5で2x²−2x−1、T03はy=3、法5で
−2x⁵−2x⁴+2x³−x²+2x−1となり、次数を保つ既約多項式である。
小さいRabin判定（x^(5^n)=x、gcd(x^5−x,g)=1）を係数演算で検査した。
Gaussの補題と次数を保つ特殊化により、元fのQQ[x,y]既約性が従う。
これは有限体標本から複素点の所属を外挿する用途ではなく、**既約性の厳密な有限証明**である。

rho=N/Dとし、曲線の接微分delta=f_y*∂x−f_x*∂yを使う。
Vandermonde、D,N,N−D,N+D、delta(rho)の分子、およびN²+D²の各剰余が非零であることを確認した。
QQ[x,y]/(f)は整域なので、これらの**積も非零**である。
従って複素数へ広げた少なくとも一つの幾何成分で、同じ非空開上で
固有値が相異なり、s(s+3)≠0、lambda*mu≠0、delta(a)≠0が同時に成り立つ。
絶対既約性を無断で仮定していない。

実際gammaのPfaffianは

    2rho*x*y*(x−1)(x−2)(y−1)(y−2)(y−x)

で、その開で非零である。
一般の曲線点ごとにt≠0のPi_tは実際のT02またはT03に入り、t=0のQaも同じ閉成分に入る。
a-imageがA1で稠密なので、Qa写像のT_i逆像が閉であることから**全有限複素a**の所属が従う。
a=0やcoverの特別な値もこの閉包に含める。

## 6. 全120labelの固定arc試験と旧T37所属

最初の固定D=diag(0,1,2,3,4)では、T01以外に得た十分条件は次の四つだった：

| 成分 | rho | a |
|---|---:|---:|
|T02|−11/28|2464/5967|
|T02|−8/25|1600/5049|
|T03|−11/29|319/810|
|T03|−8/29|1856/6993|

各値に4置換、全40×120条件を保存した。
一つだけ固有値を自由にしたdiag(0,1,2,3,d)のQQ(d)上では、T01以外に非定数rho(d)の解は無かった。
これはこの検査法の結果であって、他成分への非所属ではない。
二固有値の曲線で得た全a所属を、この四値だけへ限定してはいけない。

旧Qa1∈T37には実際の別証拠がある。
v23の `P(b)+t*tau(b)` とMPS Example5.6の全括弧一致が保存されており、
tau(1)の四つの非零上成分はQa1と正確に一致する。
t→∞の射影極限でQa1∈T37となる。
これは今回の固定toric arcのperiod試験から得た結果ではない。
T37は内部IDであり、v23の文献同定と元certificateを参照する。

producer: `code/q_mps_toric_arc_memberships_20261103_v63.py`。
曲線係数: `code/data/q_mps_toric_arc_curves_20261103_v63.json`。
保存先: `reports/Q_MPS_TORIC_ARC_MEMBERSHIPS_20261103_v63/`。
新成分数、全局所枝、他のMPS成分への非所属、全Poisson空間の尽くしは主張しない。
