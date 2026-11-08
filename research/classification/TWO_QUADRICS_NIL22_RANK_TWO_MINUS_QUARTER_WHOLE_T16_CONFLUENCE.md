# nil22長translation rank二の(-1/4,0)全core：真T16への実衝突

2026-11-08、PRIMARY E、純手全parameter曲線。新operator、point、rank、minor、profile、
jet又はCASは実行していない。同じ有限恒等式のreadbackは別receiptとする。
旧全source／独立証明書／主稿／global台帳／Gitは不変。新成分増分は零。

**命題。** `TWO_QUADRICS_NIL22_LONG_FIVE_WHOLE_SOURCES.md` §5の全core
\((\eta,\xi)=(-1/4,0)\)、全 \(s\in\mathbb C\) のactual HQQ sourceはnamed T16に入る。
他の二つのrank二coreをこの命題に含めない。

## 1. 戻すべきwhole source

既独立五族還元の元二十係数substitutionと全affine逆を入力にする。
\(R=hx-y^2/2,S=hz-yw,U=R+w^2/8\) と置けば、全六括弧
\((xy,xz,xw,yz,yw,zw)\) は正確に
\[
 P_s=(U/2+shw/4,\ s(y^2+w^2/4),\ S/2+shy,
                    -S/2+shy,\ sh^2,\ -2U+shw).             \tag{1}
\]
h行は全零、cofactorは \(y\partial_x+h\partial_y+w\partial_z\)。
actual divisorは \(-h(U^2+S^2/4)\)。
二factor \(U\pm iS/2\)、\(i^2=-1\) は相異rank三quadric。
s非零だけに制限せず、その零補集合を同じ曲線に保持する。

全六式は既potentialからも直接得られる。
\(u=x-y^2/2,v=z-wy\) において
\(\Phi=-U^2-S^2/4\)、
\(Vu=-u/2-w^2/16-sw/4,Vv=-v/2,Vw=s\)。
元座標へ戻す \(P=\partial_y\wedge V+\operatorname{Jac}\Phi\) の六係数が(1)である。
補助の非線形Ga座標をPGL同値として使わない。

## 2. 真T16の全period正逆とEuler gauge

入力は `reports/MPS_S61_RAW_CATALOGUE_20261103_v62/T16.json` の全十entry。
source座標を \((A,B,C,D,H)=(x_0,x_1,x_2,x_3,x_4)\)、
periodを \((\rho_0,\rho_1,\rho_2)=(p_0,p_1,p_2)\)、edgeを \((\tau_0,\tau_1)=(t_0,t_1)\) とする。
`TWO_QUADRICS_CONFLUENCE_NAMED_MPS_PARENTS.md` の全六同定を用い、
\[
 p=(\rho_2-\rho_1)/2,\quad q=\rho_0+(\rho_1+\rho_2)/2,
 \quad k=-(\rho_1+\rho_2)/4
\]
とすると、Euler場から
\[
 L=(0,-pB,-kC,kD,0)
\]
を引いたH-Casimir gaugeの全六式は
\[
\begin{array}{lll}
 P_{AB}=pAB+\tau_0H^2,&P_{AC}=kAC,&P_{AD}=-kAD,\\
 P_{BC}=-kBC,&P_{BD}=kBD,&P_{CD}=qCD+\tau_1H^2.
\end{array}                                                   \tag{2}
\]
実逆は
\[
 (\rho_0,\rho_1,\rho_2,\tau_0,\tau_1)
 =(q+2k,-p-2k,p-2k,\tau_0,\tau_1).                           \tag{3}
\]
raw period Pfaffian openは \(pq\ne0\)。
sourceの閉包がnamed T16であることは既受理MPS catalogueへの依存である。

## 3. 二つのparabolic frameと中間の全六多項式

\(i^2=-1\)、eをcurve parameter、sを全source parameterとする。
\[
 p=q=\tau_0=\tau_1=-e,\qquad k=is e^2.                       \tag{4}
\]
従って真raw bindingは
\[
 (\rho_0,\rho_1,\rho_2,\tau_0,\tau_1)
 =(-e+2is e^2,e-2is e^2,-e-2is e^2,-e,-e).
\]
period openは \(e^2\ne0\)、s零でも変わらない。

中間target座標を \((h,a,Y,b,W)\) とし、sourceを
\[
 A=h+eY+e^2a,\quad B=-h+eY-e^2a,\quad
 C=h+eW+e^2b,\quad D=-h+eW-e^2b,\quad H=h                 \tag{5}
\]
と置く。このsource-from-target frameのdeterminantは \(4e^6\)、全逆は
\[
 h=H,\quad a=(A-B-2H)/(2e^2),\quad Y=(A+B)/(2e),\quad
 b=(C-D-2H)/(2e^2),\quad W=(C+D)/(2e).
\]
\(R_a=ha-Y^2/2,R_b=hb-W^2/2,H_a=h+e^2a,H_b=h+e^2b\) と書く。
(2)のactual逆pushforwardは正確に
\[
\begin{array}{lll}
 P_{aY}=R_a+e^2a^2/2,&P_{ab}=isYW,&P_{aW}=isYH_b,\\
 P_{Yb}=isH_aW,&P_{YW}=isH_aH_b,&P_{bW}=R_b+e^2b^2/2.
\end{array}                                                   \tag{6}
\]
例えば \(AB+h^2=e^2(Y^2-2ha-e^2a^2)\)、
\(CD+h^2=e^2(W^2-2hb-e^2b^2)\) であり、全inverse-frame poleが相殺する。

## 4. 固定actual mapと最終全六式

さらにdeterminant \(-1\) の固定map
\[
 x=(a+b)/2,\quad y=(Y+W)/2,\quad
 z=i(a-b),\quad w=i(Y-W),\quad h=h                         \tag{7}
\]
を掛ける。全実逆は
\[
 a=x-iz/2,\quad b=x+iz/2,\quad
 Y=y-iw/2,\quad W=y+iw/2.
\]
全source-from-final frameのdeterminantは \(-4e^6\)、
actual source-to-finalのdeterminantは \(-1/(4e^6)\)。両frame逆を実際に戻せる。

(7)の右辺のa,b,Y,Wを以後このfinal座標のlinear式として使い、
\(A_e=R_a+e^2a^2/2,B_e=R_b+e^2b^2/2,K=is\) と略記する。
最終全六多項式は
\[
\begin{aligned}
 Q_{xy}&=(A_e+B_e+KYH_b-KH_aW)/4,\\
 Q_{xz}&=sYW,\\
 Q_{xw}&=i(A_e-B_e-KYH_b-KH_aW)/2,\\
 Q_{yz}&=i(-A_e+B_e-KH_aW-KYH_b)/2,\\
 Q_{yw}&=sH_aH_b,\\
 Q_{zw}&=-A_e-B_e+KYH_b-KH_aW.
\end{aligned}                                                \tag{8}
\]
h行は全零。全spatial degree二、全e冪非負、最大四である。
全e非零では、真raw T16から(2)のEuler gaugeと二つのactual mapを合成したもの。
rawとの差は同じmapでpushした \(E\wedge L\) のみであり、rawとgaugeを同一視しない。

中心では \(A_e+B_e=2U,A_e-B_e=-iS\)、
\(Y-W=-iw,Y+W=2y,K=is\) だから、(8)の全六式は厳密に(1)となる。
特に全磁気sは \(Q_{yw}=sh^2\) に残り、一点又はs零sliceから全源を外挿していない。

## 5. 全actual Pfaffianとpointwise結論

\(f_1=Y^2-2ha-e^2a^2,f_2=W^2-2hb-e^2b^2\) とすると、
sourceのFとGはそれぞれ \(-e^3f_1,-e^3f_2\)。
(2)の全PfaffianはFGであり、交差k項は正確に相殺する。
従って全actual divisorは
\[
                       D_e=-h f_1f_2/4.                     \tag{9}
\]
中心は \(-h(U^2+S^2/4)\ne0\) である。
(8)は全eの実多項式curveで、e非零は真named T16。
閉性から全sの中心もT16へ入り、既full20 sourceの全affine逆を戻して
**rank二(-1/4,0)whole sourceはT16**。

型七active-u／active-v、他のrank二二core、他の因子型は独立に扱う。
本稿保存時には新しい独立算術receiptを主張しない。
後続のtiny QQ(i)-Laurent readbackは(2)–(9)の同じ全e/s曲線だけを確認し、
新operator、point、rank又はprofileを追加しない。
