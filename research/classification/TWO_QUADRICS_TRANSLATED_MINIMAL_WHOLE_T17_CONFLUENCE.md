# H＋二既約quadric・最小冪零長translation：全源の真T17実衝突

2026-11-08、PRIMARY `E_reduced_nonnormal`。純手で導いた同一の全parameter曲線である。
本稿保存時に新operator、point、rank、minor、profile、jet又はCASは実行していない。
後続の有限Laurent readbackは別receiptとする。新成分増分は零。
旧producer／immutable証明書／主稿／global台帳／Gitは変更しない。

**命題。** `TWO_QUADRICS_PURE_NILPOTENT_SUPPORT.md` の型十、
cofactor \(Z=y\partial_x+\partial_y\) の全 \(H+q_1+q_2\) actual sourceはnamed T17に入る。
\(q_1,q_2\) は相異なる既約quadricであり、以下では全磁気translation \(\kappa\) を零にも非零にも制限しない。

## 1. 全sourceと正確な因子条件

`HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_NORMAL_FORMS.md` §§1–2の全二十五係数potential表示を用いる。
その元全九十係数との可逆affine map、全Jacobi、完成平方 \(R=r+L/\delta\) の実逆を保持する。
HQQ因子条件による全分岐は `TWO_QUADRICS_PURE_NILPOTENT_SUPPORT.md` §6と
`TWO_QUADRICS_PURE_NIL_1_4_10_WHOLE_REDUCTION.md` §3であり、既約quarticの帰属を流用しない。
\(h=1\) の補助 \(r=x-y^2/2\) で \(\Phi=\delta r^2/2+az^2\)、
\(V=-r\partial_r/2-z\partial_z/2+\kappa\partial_w\) となる。
\(\kappa\) は任意、\(\delta a\ne0\) が正確な対象条件である。

斉次座標順を \((h,x,y,z,w)\)、\(R=hx-y^2/2\) とすると、全六括弧の順序
\((xy,xz,xw,yz,yw,zw)\) において
\[
 P(\delta,a,\kappa)=
 (R/2,-yz/2,-2ahz+\kappa hy,-hz/2,\kappa h^2,\delta R).       \tag{1}
\]
h行は全零。actual divisorは
\[
 D=h\left(\delta R^2/2+ah^2z^2\right).
\]
\(\delta=0\) は最高quartic零、\(a=0\) は平方で対象外。
\(\delta a\ne0\) では二因子 \(R\pm c hz\)、\(c^2=-2a/\delta\) は
相異なるrank三quadricであり、hで割れない。
元potentialの完成平方は \(x\mapsto x+L(z,w)/\delta\) という実affine mapなので、
補助の非線形Ga座標をPGL同値へ取り違えていない。

## 2. 真completed T17の全period couplingとEuler gauge

sourceは `reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T17.json` の全十entryである。
原座標を \((A,B,C,D,H)=(x_0,x_1,x_2,x_3,x_4)\)、periodを
\((\rho_0,\rho_1,\rho_2)=(p_0,p_1,p_2)\)、smoothingを \((\tau_0,\tau_1)=(t_0,t_1)\) とする。
特に raw \(P_{CD}\) の \(\tau_0\tau_1H^2/\rho_1\) を除かない。
Euler場を \(E\) とし、
\[
 L=(0,-pB,-\alpha C,-\beta D,0)
\]
を引くと \(P_T-E\wedge L\) はH-Casimirとなり、全六括弧は
\[
 P_{AB}=F,\quad P_{AC}=\alpha AC,\quad P_{AD}=\beta AD,\quad
 P_{BC}=-\alpha BC,\quad P_{BD}=-\beta BD,\quad P_{CD}=G,
\]
\[
 F=pAB+uH^2,\qquad G=qCD+vF,\qquad \alpha+\beta=p.           \tag{2}
\]
これは \(vAB\) ではなく **\(vF\) 全体**を持つcoupled sourceである。
`TWO_QUADRICS_CONFLUENCE_NAMED_MPS_PARENTS.md` の明示正逆は
\[
 p=\rho_1,\quad \alpha=(\rho_1-\rho_2)/2,\quad
 \beta=(\rho_1+\rho_2)/2,\quad q=\rho_0+\rho_2,
 \quad u=\tau_0,\quad v=\tau_1/\rho_1,
\]
\[
 (\rho_0,\rho_1,\rho_2,\tau_0,\tau_1)
 =(q-\beta+\alpha,p,\beta-\alpha,u,pv).                       \tag{3}
\]
raw period Pfaffian openは \(\rho_1(\rho_0+\rho_2)=pq\ne0\)。
全completed sourceのnamed T17認識は既受理catalogueへの依存であり、
本稿の新部分は以下のactual全tensor同定である。

## 3. 全parameter binding、frame及び実逆

\(e\) を曲線parameter、\(\delta\rho\ne0\)、\(\kappa\) 任意とする。
まずshift前のtarget座標で、source座標を
\[
 H=h,\quad A=h+ey+e^2x,\quad B=h-ey+e^2x,\quad
 C=h+e(w+2\rho y),\quad D=z                              \tag{4}
\]
と置く。targetからsourceへのlinear frame \(J_e\) のdeterminantは \(2e^4\)。
その全実逆は
\[
 h=H,\quad x=(A+B-2H)/(2e^2),\quad y=(A-B)/(2e),\quad
 z=D,\quad w=(C-H-\rho(A-B))/e.                            \tag{5}
\]
source順は \((A,B,C,D,H)\)、target順は \((h,x,y,z,w)\) であり、
(4)の行をその順に並べたものが \(J_e\)、(5)が \(J_e^{-1}\) である。

真helperの全parameterを
\[
 p=-e/2,\quad u=e/2,\quad q=-2\rho e,\quad v=\delta/e^2,
 \quad \alpha=\kappa e^2,\quad \beta=-e/2-\kappa e^2          \tag{6}
\]
とする。従って真raw completed bindingは
\[
 (\rho_0,\rho_1,\rho_2,\tau_0,\tau_1)
 =(-2\rho e+e/2+2\kappa e^2,-e/2,-e/2-2\kappa e^2,e/2,-\delta/(2e)). \tag{7}
\]
\(\rho_1(\rho_0+\rho_2)=\rho e^2\ne0\) なので、任意の \(e\ne0\) で真period openに入る。
\(\kappa=0\) でもこの単元は変わらない。

\[
 T_e=h+e^2x,\quad R_e=hx-y^2/2+e^2x^2/2,\quad
 C_e=h+e(w+2\rho y),\quad L_e=1/2+\kappa e
\]
と書くと \(AB=h^2+2e^2R_e\) なので、**全period poleが相殺し**
\[
 F=-e^3R_e,\qquad G=-e(\delta R_e+2\rho C_ez).              \tag{8}
\]
(2)の全六式を(5)でactual pushforwardした結果は正確に
\[
 (R_e/2,-L_eyz,\kappa yC_e-\rho R_e,-L_eT_ez,
     \kappa T_eC_e,\delta R_e+2\rho(C_e-L_eT_e)z).             \tag{9}
\]
全spatial degreeは二、全e冪は非負で最大三。
単なるPfaffian一致ではなく、全十entryのh行零及び全六entryの恒等式である。

## 4. 中心をcanonical sourceへ戻すconstant x-shift

さらに実map \(x_{\rm new}=x_{\rm old}+\rho z/\delta\) を掛ける。
determinant一、逆は同じ係数を引く。全 \(\delta\ne0\) で正則である。
最終座標で
\[
 X=x-\rho z/\delta,\quad R_e=hX-y^2/2+e^2X^2/2,
 \quad T_e=h+e^2X,\quad C_e=h+e(w+2\rho y),\quad L_e=1/2+\kappa e
\]
と定義し直すと、最終全六式は
\[
\begin{aligned}
 Q_{xy}&=R_e/2+(\rho/\delta)L_eT_ez,&Q_{xz}&=-L_eyz,\\
 Q_{xw}&=\kappa yC_e+(2\rho^2/\delta)(C_e-L_eT_e)z,
 &Q_{yz}&=-L_eT_ez,\\
 Q_{yw}&=\kappa T_eC_e,&Q_{zw}&=\delta R_e+2\rho(C_e-L_eT_e)z.
\end{aligned}                                                  \tag{10}
\]
h行は全零。generic eで真T17からの全actual mapは
\(J_e^{-1}\) にこのconstant shiftを合成したもの、determinantは \(1/(2e^4)\)。
その逆は (4) の x を X に置いた全五式であり、両積はidentityである。
raw sourceとの差は同じmapでpushした \(E\wedge L\) だけ。
このEuler項を省いてraw completionそのものを(10)と主張していない。

h-Casimir gaugeのaffine Pfaffian、従ってactual divisorは
\[
 \operatorname{Pf}_{4}(Q_e)=R_e(\delta R_e/2+\rho C_ez),\qquad
                 D_e=hR_e(\delta R_e/2+\rho C_ez).             \tag{11}
\]
実determinant変換則でも (2) の \(FG\) と \(\det J_e=2e^4\) から同じ式が出る。
これは全126 quintic係数を持つ式で、最高部だけではない。

中心e零で \(R_e=R-\rho hz/\delta\)、\(T_e=C_e=h,L_e=1/2\)。
(10)は厳密に(1)、ただし
\[
                       a=-\rho^2/(2\delta)                    \tag{12}
\]
となる。特に \(Q_{xw}=\kappa hy+\rho^2h z/\delta\)、
\(Q_{zw}=\delta R\)。actual Pfaffianは
\(\delta R^2/2-\rho^2h^2z^2/(2\delta)\) である。
cofactorも \(Z=y\partial_x+h\partial_y\) に戻る。

## 5. 全unit／零補集合とpointwise結論

任意の対象 \(\delta a\ne0\)、任意 \(\kappa\) に対し、
\(\rho^2=-2\delta a\) の二つの非零複素平方根を選ぶ。
\(\kappa\) は反転せず全零補集合を同じ式に残し、対象の \(\delta,a\) 以外のunitを仮定しない。
(10)はe零で非零の実多項式曲線であり、e非零の全点は真completed T17。
よってnamed T17の閉性から(1)もT17に入る。
元二十五係数への全affine完成平方／固有座標／scalar／PGLの明示逆を戻せば、
**型十の全HQQ actual sourceがT17に入る**。

既 \(\kappa=0\) の別MPS曲線との両立は許され、その部分が複数閉成分の交点に入っても矛盾しない。
旧partial noteの「\(\kappa\ne0\) の親未認識」は作成時の状態として変更しない。
型七の残四軌道、他の因子型、全分類は本稿の結論に含めない。
有限readbackは(4)–(11)の同一全parameter恒等式だけを対象とし、
新anchor、operator、rank、profile、higher jet又は別curveを追加しない。
