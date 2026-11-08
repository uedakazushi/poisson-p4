# nil22長translation active-v全源：真T18へのprincipal-nil実衝突

2026-11-08、PRIMARY E。純手で同一の全parameter曲線を構成する。
新operator、point、rank、minor、profile、jet又はCASは実行していない。
後続の同じ有限恒等式のreadbackは別receiptとする。新component増分は零。
旧source／独立五族証明／主稿／global台帳／Gitは変更しない。

**命題。** `TWO_QUADRICS_NIL22_LONG_FIVE_WHOLE_SOURCES.md` §4のwhole active-v族、
全 \(c\ne0,B\) はnamed T18に入る。
rank一active-uの既T17曲線、rank二三core及び他因子型とは区別する。

## 1. 対象whole sourceと全元係数への逆

元全二十係数potentialと実normalizerは既五族純証明とその独立査読を入力にする。
cofactorは \(Z=y\partial_x+h\partial_y+w\partial_z\)、\(R=hx-y^2/2\)、\(S=hz-yw\)。
全六括弧の順序 \((xy,xz,xw,yz,yw,zw)\) で、対象は正確に
\[
 P(c,B)=\bigl(Bhw+cw^2+R/4,\ xw/4-3yz/8+Bw^2,
 (hz-3yw)/8,\ yw/4-hz/2,\ -hw/4,\ -w^2/4\bigr)           \tag{1}
\]
を意味する。h行は全零。
actual divisorは \(-h(S^2+4cw^4)/16\)。
\(c\ne0\) が相異二rank四quadricの正確な条件、\(c=0\) は平方として対象外である。
B零は除かず、全元lower／magnetic係数は既actual affine逆で戻せる。

## 2. 真completed T18とh-Casimir gauge

入力は `reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T18.json` の全十entry。
source座標を \((A,B_s,C,D,H)=(x_0,x_1,x_2,x_3,x_4)\)、
periodを \((p_0,p_1)=(r,s)\)、smoothingを \((t_0,t_1)=(t,u)\) とする。
特にraw \(P_{AC}=-2tuH^2/(2r-s)\) の完成項を省かない。
Euler場をEとし、
\[
 L=(0,(r-s)B_s,(r-s/2)C,-sD/2,0)
\]
を引くとH-Casimir gaugeになる。
\(d=r-s,\nu=r-s/2\) と書くと、その全六式は
\[
\begin{array}{lll}
 P_{AB_s}=-dAB_s+tDH,&P_{AC}=-\nu AC-tuH^2/\nu,&
 P_{AD}=(\nu-d)AD,\\
 P_{B_sC}=(\nu-d)B_sC,&P_{B_sD}=-(\nu-d)B_sD,&
 P_{CD}=dCD+uB_sH.
\end{array}                                                   \tag{2}
\]
periodの正逆は \(r=2\nu-d,s=2\nu-2d\)。
raw period Pfaffian openは \(-d\nu\ne0\)。
全named T18認識は既受理catalogueへの依存であり、
本稿では因子型又はcofactor固有値だけからmembershipを推測しない。

## 3. 全parameter bindingとframeの正逆

\(e\) をcurve parameter、\(\alpha\ne0\)、\(\beta\) 任意とし、
\[
 k=1/4+\beta e,\quad d=e/2,\quad\nu=ek,\quad
 (r,s,t,u)=(2\beta e^2,-e/2+2\beta e^2,\nu,-\nu).           \tag{3}
\]
\(tu=-\nu^2\) なので、真raw correctionは正確に \(\nu H^2\) へ相殺し、
(2)の \(P_{AC}=\nu(H^2-AC)\) となる。
\(e k\ne0\) の非空稠密開では真sourceのperiod open \(-e^2k/2\ne0\) に入る。
\(\beta=0\) もこの同じ開と式に含まれる。

source座標をtarget \((h,x,y,z,w)\) から
\[
\begin{aligned}
 H&=h,&A&=h+ey+\alpha ew+e^2x,&B_s&=w-ez/2,\\
 C&=h-ey+\alpha ew+e^2x,&D&=w+ez/2
\end{aligned}                                                \tag{4}
\]
と置く。source順 \((A,B_s,C,D,H)\) のframe determinantは \(-2e^4\)。
実逆は
\[
 h=H,\quad y=(A-C)/(2e),\quad z=(D-B_s)/e,\quad
 w=(B_s+D)/2,\quad
 x=(A+C-2H-\alpha e(B_s+D))/(2e^2).                        \tag{5}
\]
全e非零でactual PGLであり、raw Euler項にも同じ実mapを用いる。
\(1/e\) を持つmap自体と、以下の全e正則な曲線を混同しない。

## 4. 全六多項式と中心の全係数

\(T=h+\alpha ew+e^2x\)、\(R=hx-y^2/2\) とする。
(2)–(5)のactual pushforwardは正確に
\[
\begin{aligned}
 Q_{xy}={}&kR+2\alpha\beta hw+\alpha^2(k-1/2)w^2
 +\alpha e(3k-1)xw/2+ke^2x^2/2+\alpha ke yz/4,\\
 Q_{xz}={}&kxw+(k-1)yz/2+2\alpha\beta w^2
                  -\alpha e(k-1/2)z^2/4,\\
 Q_{xw}={}&(k-1)yw/2+k(T+h)z/4,\\
 Q_{yz}={}&((k-1)T-kh)z/2+kyw,\\
 Q_{yw}={}&((k-1)T+kh)w/2+ke^2yz/4,\\
 Q_{zw}={}&(k-1/2)(w^2-e^2z^2/4).
\end{aligned}                                                \tag{6}
\]
h行は全零。全spatial次数は二、全e冪は非負、最大三である。
例えば \(Q_{xy}\) の一見した \(\alpha(2k-1/2)hw/e\) は
\(2\alpha\beta hw\) へ厳密に相殺する。
この項がBを全て動かし、固定period比のcurveだけではB零しか出ない点が本質である。

中心では \(k=1/4,T=h\) なので、(6)の全六係数は厳密に(1)で
\[
                         c=-\alpha^2/4,\qquad B=2\alpha\beta. \tag{7}
\]
\(\alpha\) の二つの非零複素平方根で任意c非零を覆い、
\(\beta=B/(2\alpha)\) が全Bを覆う。B又はβを反転する必要は無い。

## 5. Actual Pfaffianと全sourceへの閉性

(2)の全affine Pfaffianは、交差項を含め
\[
 -d\nu\,(AB_s-tDH/\nu)(CD+uB_sH/\nu).
\]
(3)では二因子 \(Q_1=AB_s-DH\)、\(Q_2=CD-B_sH\) となり、
\[
\begin{aligned}
 Q_1/e={}&-S+\alpha w^2
       +e\{xw-(\alpha w+y)z/2\}-e^2xz/2,\\
 Q_2/e={}& S+\alpha w^2
       +e\{xw+(\alpha w-y)z/2\}+e^2xz/2.
\end{aligned}                                                \tag{8}
\]
actual mapのdeterminant \(-1/(2e^4)\) を掛けると、全actual divisorは
\[
                     D_e=h\,k(Q_1/e)(Q_2/e)/4.                \tag{9}
\]
中心は \(-h(S^2-\alpha^2w^4)/16\)、\(\alpha\ne0\) で非零であり、
対象の二既約quadric積そのもの。

(6)は全eで正則、generic \(ek\ne0\) では真completed T18の全actual像。
従って全Jacobiは多項式恒等式として中心にも成り立ち、閉性から中心もT18に入る。
(7)と元full20 sourceの全actual逆を戻せば、**active-v whole族、全c非零／全BはT18**。
rank二三coreのs非零軌道とその零境界は本稿で解決したとはしない。

有限readbackは、(2)–(9)の同じ全e/α/β曲線だけを確認する。
raw nonmonomial denominator \(2r-s=2\nu\) は
\((2r-s)P_{AC}=-2tuH^2\) をwhole係数で掛け戻して検証でき、
根k零で親sourceを定義したと仮定しない。
別anchor、operator、profile又はhigher jetを追加する必要は無い。
