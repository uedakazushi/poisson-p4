# nil22長translation rank二(0,0)全core：moving-H真T18実衝突

2026-11-08、PRIMARY E。全磁気sを持つ同一の純手曲線を構成する。
新operator、point、rank、minor、profile、GB、jet又はCASは実行していない。
後続の有限readbackは別receipt。旧全source／主稿／global台帳／Gitは不変。
新component増分は零。

**命題。** `TWO_QUADRICS_NIL22_LONG_FIVE_WHOLE_SOURCES.md` §5の
全 \((\eta,\xi)=(0,0)\)、全 \(s\in\mathbb C\) actual HQQ sourceはnamed T18に入る。
親のhyperplaneは曲線上で動く。hをgeneric曲線のCasimirと仮定しない。

## 1. 対象whole source

既独立五族還元の全二十係数substitution、元全potential及びactual affine逆を入力とする。
\(R=hx-y^2/2,S=hz-yw\) とすると、全六括弧の順序
\((xy,xz,xw,yz,yw,zw)\) で対象は
\[
 P_s=(R/2,\ xw/2-yz/4+sy^2,\ S/4+shy,
                -S/2+shy,\ sh^2,\ -2R+shw).                \tag{1}
\]
h行零、cofactor \(y\partial_x+h\partial_y+w\partial_z\)、
actual divisor \(-h(R^2+S^2/8)\)。二factorは相異rank四quadric。
s零を含め全sourceを扱う。他のrank二二coreの帰属を外挿しない。

## 2. 真completed T18、全period及び動くhyperplane

真sourceは `reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T18.json` の全十entry、
source座標は \((A,B,C,D,H)=(x_0,x_1,x_2,x_3,x_4)\)。
periodを \((p_0,p_1)=(r_p,s_p)\)、edgeを \((t_0,t_1)=(t,u)\) とする。
\(r_0^2=-1/2\)、全sに対して \(\sigma=2r_0s\) と置き、
\[
 r_p=er_0,\quad s_p=e^2\sigma,\quad
 \nu=r_p-s_p/2=e(r_0-e\sigma/2),\quad t=\nu,\quad u=-\nu.  \tag{2}
\]
真raw correction \(-2tuH^2/(2r_p-s_p)\) は正確に \(\nu H^2\)。
raw period openは
\[
 -(r_p-s_p)(r_p-s_p/2)
     =-e^2(r_0-e\sigma)(r_0-e\sigma/2)\ne0.                \tag{3}
\]
この非空稠密開だけで親sourceを使う。二つの因子が零になる根へ元分母を代入しない。
全eの曲線は下の多項式表示で延長する。

中間target座標を \((h,a,Y,c,W)\) として
\[
 A=h+e(W-Y)+e^2a,\quad B=h+eY,\quad
 C=h+e(Y+W)+e^2c,\quad D=h,\quad H=h+eW.                 \tag{4}
\]
source-from-target frameのdeterminantは \(-e^6\)。全逆は
\[
 h=D,\quad a=(A+B-H-D)/e^2,\quad Y=(B-D)/e,
 \quad c=(C-B-H+D)/e^2,\quad W=(H-D)/e.                    \tag{5}
\]
親のactual factorはHであり、(4)ではH=h+eW、中心でhへ動く。
Hをgeneric全曲線のhと同一視しない。

## 3. 全actual map後のEuler gaugeと全多項式

(5)でraw全tensorをactual pushforwardしたものをPrawとする。
中間Euler場Eに対して次の全linear場を引く。
\[
\begin{aligned}
 L_h&=0,\quad L_a=r_0h/e-\sigma h/2+r_0Y-e\sigma W/2,\\
 L_Y&=r_0(h+eY),\quad
 L_c=e\sigma Y/2+r_0W+e^2\sigma c,\quad
 L_W=e\sigma(h+eW)/2.
\end{aligned}                                               \tag{6}
\]
L自体のpoleを勝手に無視せず、\(Q=Praw-E\wedge L\) の全係数で相殺する。
\(H_e=h+eW\) と書けば、唯一の非零h行は
\[
                         Q_{hc}=e(r_0-e\sigma/2)YW.          \tag{7}
\]
順序 \((aY,ac,aW,Yc,YW,cW)\) の残全六式は
\[
\begin{aligned}
 Q_{aY}={}&-r_0(ha+YW-Y^2)+\sigma hY/2
                              +e(\sigma YW/2-r_0aY),\\
 Q_{ac}={}&r_0\{Yc-W(c+a)\}+\sigma Y^2
                      +e\sigma Y(c-a/2)-e^2\sigma ac,\\
 Q_{aW}={}&\sigma\{H_eY-hW/2-e(aH_e+W^2)/2\},\\
 Q_{Yc}={}&r_0(hc-YW)+\sigma hY/2
                      +e\{\sigma Y(Y+W)/2+r_0cY\},\\
 Q_{YW}={}&\sigma H_e(h+eY)/2,\\
 Q_{cW}={}&r_0(YW-hc)+\sigma hW/2-er_0cW
              +e\sigma(-YW/2+W^2/2+hc)+e^2\sigma cW.
\end{aligned}                                               \tag{8}
\]
全spatial degree二、全e冪非負、最大二である。
(7)を落としてh-Casimir familyと主張してはならない。
全generic eではprojective Jacobi \(E\wedge[Q,Q]=0\) が成り立つ。
strict Jacobiをgeneric eで要求する必要はなく、中心ではh行零から全strict Jacobiが従う。

## 4. 固定実frameとcanonical中心の全六係数

次の固定actual mapを掛ける。
\[
 x=(a+c)/(8r_0^2),\quad y=-Y/(2r_0),\quad
 z=(c-a)/(2r_0),\quad w=Y-2W,\quad h=h.                   \tag{9}
\]
\(r_0^2=-1/2\) を使うと実逆は
\[
 a=-2x-r_0z,\quad c=-2x+r_0z,\quad
 Y=-2r_0y,\quad W=-r_0y-w/2.
\]
固定determinantは \(1/(8r_0^4)=1/2\)、
全source-from-final frameは \(-2e^6\)、actual逆は \(-1/(2e^6)\)。
全raw Euler項にも同じ実mapとinverseを用いる。

(8)の中心を \(\sigma=2r_0s\) とともに(9)でpushすると、全六係数は正確に(1)。
例えば \(q_1=ha+YW-Y^2,q_2=hc-YW\) について
\[
 q_1=4r_0^2R-r_0S,\qquad q_2=4r_0^2R+r_0S.
\]
\(Q_{xy}=R/2,Q_{xz}=xw/2-yz/4+sy^2,Q_{yw}=sh^2\) であり、
残三式も(1)と一致する。sを反転又は一点へ固定していない。

## 5. 全actual Pfaffianとpointwise結論

(2)の真T18の二factorは \(AB-DH\)、\(CD-BH\) であり、(4)では厳密に
\[
 AB-DH=e^2f_1,\quad f_1=ha+YW-Y^2+eaY,
 \qquad CD-BH=e^2f_2,\quad f_2=hc-YW.
\]
raw actual divisorは \(-(r_p-s_p)\nu H(AB-DH)(CD-BH)\)。
(5)のactual determinant \(-e^{-6}\) を掛けると、中間の全actual divisorは
\[
           H_e(r_0-e\sigma)(r_0-e\sigma/2)f_1f_2.            \tag{10}
\]
Euler gaugeはactual divisorを変えない。さらに(9)のdeterminant1/2を掛け、
中心で \(r_0^2=-1/2\) を使えば正確に \(-h(R^2+S^2/8)\) となる。
これは全126係数の恒等式で、最高部又はselected factorだけの比較ではない。

(7)–(8)と固定mapは全eで多項式、中心は非零でgeneric eは真completed T18。
したがって閉性から全sの中心もnamed T18に入り、元full20の全actual逆を戻して
**rank二(0,0)whole sourceはT18**。

\(r_0\) の二代数選択は定数frameの選択であり、近似複素数は使わない。
s零を含め全parameterを同じ式に保持する。
他のrank二core、他因子型又は全分類の完了をこの証明から推論しない。
後続のtiny QQ(\(\sqrt{-2}\))-Laurent readbackは同じ全e/s曲線だけを確認し、
新operator、point、profile又はhigher jetを追加しない。
