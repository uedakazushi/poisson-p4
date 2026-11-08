# nil22長translation最後のrank二core：quadric-pencilの真Ga²による全T19実衝突

2026-11-08、PRIMARY E。新しい構造的whole-source曲線を純手で構成する。
既D102／M47のlocal routeは停止し、追加obstruction／jet／point／profileは無い。
本稿保存時に新operator又はCASは実行していない。
後続の同一有限恒等式readbackは別receipt。新成分増分は零。
旧全source／不成立frame診断／主稿／global台帳／Gitは変更しない。

**命題。** `TWO_QUADRICS_NIL22_LONG_FIVE_WHOLE_SOURCES.md` §5の最後のwhole core
\(\eta=1/12,\xi^2=-1/432\)、全sはnamed T19に入る。
xi二符号は既actual r=−1の逆で交換されるため、以下では \(j^2=-3,\xi=-j/36\) とする。
s非零だけでなく全零補集合を同じcurveで保持する。

## 1. 真completed T19と二quadric

source座標を \((a,b,c,d,H)=(x_0,x_1,x_2,x_3,x_4)\)、
periodを \((p_0,p_1)=(r,v)\)、edgeを \((t_0,t_1)=(t,u)\) とする。
真入力は `reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T19.json` の全十entry、
特に \(P_{ac}=2tu bH/(r-v)\) のcoupled correctionを保持する。
\[
             t=(r-v)/2,\qquad u=r-v.                         \tag{2}
\]
よってraw completionは正確に \((r-v)bH\)、period openは \(-r(r-v)\ne0\)。
raw sourceからEuler場に
\[
 L=(0,-rb,(-3r+v)c/2,(-r-v)d/2,0)
\]
を引いたH-Casimir gaugeの全六式は
\[
\begin{array}{ll}
 P_{ab}=rab+(r-v)dH/2,&P_{ac}=(3r-v)ac/2+(r-v)bH,\\
 P_{ad}=(r+v)ad/2,&P_{bc}=(r+v)bc/2,\\
 P_{bd}=-(r+v)bd/2,&P_{cd}=-2rcd+(r-v)b^2.
\end{array}                                                   \tag{3}
\]
そのactual divisorは
\[
 D_T=-r(r-v)H q_3q_4,\quad q_3=cd-b^2,\quad q_4=ab+dH.      \tag{4}
\]
q3はrank三、q4はrank四。この順序は以後固定する。

## 2. 真のcommuting Ga² actual automorphismと全実逆

二つのlinear nilpotent場
\[
 X=d\partial_b+2b\partial_c-a\partial_H,\qquad
 Y=-d\partial_a+b\partial_H
\]
は可換であり、q3,q4をそれぞれ保存する。
全parameter \(\lambda,\mu\) の実flowは
\[
\begin{aligned}
 a'&=a-\mu d,& b'&=b+\lambda d,\\
 c'&=c+2\lambda b+\lambda^2d,&d'&=d,\\
 H'&=H-\lambda a+\mu b+\lambda\mu d.
\end{aligned}                                                \tag{5}
\]
determinant一、逆は \((\lambda,\mu)\mapsto(-\lambda,-\mu)\)。
全二quadricはliteralに不変であり、一つのpoint又はsliceの対称性を仮定しない。
以下 \(\lambda=\mu=1/e\)。実divisorの変換で使うのは**inverse pullback**
\[
 H_{\rm old}=H+(a-b)/e+d/e^2,                              \tag{6}
\]
である。forward mapのH'と(6)の符号を取り違えない。

## 3. d-chartの補助座標と全六多項式

h=dを最終chart hyperplaneとするが、generic eでは親のhyperplaneは(6)で動く。
d一のchartで計算専用の補助座標を
\[
 U=b,\quad V=-a,\quad I=c-b^2,\quad J=H+ab,
 \qquad (a,b,c,H)=(-V,U,I+U^2,J+UV)                        \tag{7}
\]
とする。これは多項式自己同型だが、PGL変換としては一切使わない。
actual global mapは(5)というlinear frameだけである。

(3)のprojective bivectorをd-chartに制限し、(7)へ変換すると、全括弧は厳密に
\[
\begin{array}{ll}
 \{U,V\}=(r+v)UV/2+(r-v)J/2,&\{U,I\}=-(r-v)UI,\\
 \{V,I\}=-(r+v)VI,&\{U,J\}=vUJ,\\
 \{V,J\}=-vVJ,&\{I,J\}=2rIJ.
\end{array}                                                \tag{8}
\]
(5)はU,Vの実translation、I,J不変である。
よってactual pushforwardは(8)のU,VをU−1/e,V−1/eに置くだけ。

全source parameter sに対し \(j^2=-3,\sigma=4js\) と置き、
\[
              r=e,\quad v=-e+e^2\sigma.                    \tag{9}
\]
全period openは \(-e^2(2-e\sigma)\ne0\)、generic eで非空である。
このbinding後、pushされた全六括弧は正確に次の多項式となる。
\[
\begin{aligned}
 \{U,V\}_e&=\sigma/2-e\sigma(U+V)/2
                       +e^2\sigma UV/2+eJ-e^2\sigma J/2,\\
 \{U,I\}_e&=(2-e\sigma)I-(2e-e^2\sigma)UI,\\
 \{V,I\}_e&=e\sigma I-e^2\sigma VI,\\
 \{U,J\}_e&=(1-e\sigma)J-eUJ+e^2\sigma UJ,\\
 \{V,J\}_e&=(-1+e\sigma)J+eVJ-e^2\sigma VJ,\\
 \{I,J\}_e&=2eIJ.
\end{aligned}                                                \tag{10}
\]
全e負冪が消え、最大e次数二。
(7)のinverseはeに依存しないpolynomialなので、元chartの全bivector係数もe多項式である。

P4のglobal quadratic projective bivector空間のd-chartへの制限は固定のinjective linear map。
generic eのclassはその像に属し、像はlinear閉空間だから、(10)の全e多項式classもその像に属する。
従って固定linear sectionを使えば、元の全150係数における多項式curveが得られる。
ここでsourceのcanonical divergence零lift
\[
              \overline P=P_{\rm raw}-E\wedge\operatorname{div}P_{\rm raw}/5
\]
を採用すれば、GL-equivarianceにより(5)でactual pushしたものがそのsectionそのもの。
これが全150係数のe負冪無しを証明する。genericにd-Casimirを仮定したわけではない。
後続readbackではこのcanonical lift、raw Euler差、両実inverseを全係数で確認する。

## 4. 中心の元座標六式と固定canonical frame

(10)の中心を(7)で戻すと、d一の順序(a,b,c,H)で
\[
\begin{array}{ll}
 P_{ab}=\sigma/2,&P_{ac}=\sigma b,\\
 P_{aH}=H+ab-\sigma a/2,&P_{bc}=2(c-b^2),\\
 P_{bH}=H+ab+\sigma b/2,&P_{cH}=2ac+2bH+\sigma b^2.
\end{array}                                                   \tag{11}
\]
斉次化では定数をd二、一次をd倍する。
この中心はd-Casimir quadratic liftを持つが、generic eもそうだとはしない。

次の固定actual frameを掛ける。
\[
 h=d,\quad x=(c-3H)/36,\quad y=-(a+b)/4,
 \quad z=-j(c+3H)/18,\quad w=j(b-3a)/6.                     \tag{12}
\]
source順(a,b,c,d,H)から最終順(h,x,y,z,w)へのdeterminantは \(-1/216\)。
h=dを先に並べ替えた順(h,a,b,c,H)からなら \(1/216\)。
全実逆は
\[
 a=-y+jw/2,\quad b=-3y-jw/2,\quad
 c=18x+3jz,\quad d=h,\quad H=-6x+jz.                       \tag{13}
\]
\(R=hx-y^2/2,S=hz-yw\)、\(\eta=1/12,\xi=-j/36,\kappa_0=-1/6\)、
\(U_0=R-w^2/24,W_0=S-jw^2/6\) とすれば、中心の全六係数は
\[
\begin{aligned}
 P_{xy}&=U_0/2-shw/12,\\
 P_{xz}&=2xw/3+jzw/18-yz/3+s(y^2-w^2/12+jyw/3),\\
 P_{xw}&=W_0/6+shy,\\
 P_{yz}&=-W_0/2+shy+js hw/3,\\
 P_{yw}&=sh^2,\\
 P_{zw}&=-2U_0+shw.
\end{aligned}                                                \tag{14}
\]
h行零、cofactorは \(N=y\partial_x+h\partial_y+w\partial_z\)。
これは元full20 sourceの当該core/s任意の全六式そのもの。
canonical strict liftの中心は(14)から \(E\wedge N/5\) を引いたものであり、
そのEuler差を無視して二つのliftをliteral同一とは呼ばない。

## 5. 全actual divisor、全unit零補集合とnamed結論

q3,q4は(5)で厳密に不変だから、(4),(6),(9),(12)から全actual divisorは
\[
 D_e=(2-e\sigma)\{d+e(a-b)+e^2H\}q_3q_4/216,              \tag{15}
\]
ただしa,b,c,d,Hは(13)の最終座標式である。
\[
 q_3=18(R+jS/6+w^2/24),\quad
 q_4=-6(R-jS/6-w^2/8)
\]
なので、中心は厳密に
\[
 -h(R+jS/6+w^2/24)(R-jS/6-w^2/8)
                =h\{-U_0^2-W_0^2/12\}\ne0.              \tag{16}
\]
全126 actual Pfaffian係数の恒等式であり、最高部だけの比較ではない。

全sを反転しない。(10)の同じ多項式curveの中心が(14)、generic eは真completed T19。
よってnamed T19の閉性から全sの中心もT19に入る。
xiの他符号、元全affine parameterは独立五族定理の全actual逆を戻せる。
**最後のrank二coreのwhole HQQ sourceはT19**。
この陽性curveは以前の二失敗frameの陰性を反転解釈したものではなく、
quadric-pairをliteralに保存する新しいGa²構造から構成した別の全source橋である。

本稿保存時には新算術receiptを主張しない。有限readbackは同じe/s/j曲線の
全150 canonical／raw Euler正逆、全126actualPf、全periodpole相殺、
(10)／(14)及び全Jacobiだけを確認し、局所obstruction／別point／profileは追加しない。
