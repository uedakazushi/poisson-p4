# H＋二既約quadric・regular nil：全 source の T20 実衝突

2026-11-08、PRIMARY E。新しい純手 actual Poisson 曲線。
本稿作成時に新 operator、point、rank、minor、profile、jet 又は CAS は実行していない。
後続の同じ有限恒等式の独立／算術 readback は別 receipt にする。
新成分増分は零。旧 frozen notes／reports／主稿／global 台帳／Git は変更しない。

**命題。** `TWO_QUADRICS_PURE_NILPOTENT_SUPPORT.md` の型一、
regular nil cofactor・translation零の全 \(H+q_1+q_2\) actual source は既知 MPS 成分 T20 に入る。
\(q_i\) は相異なる既約 quadric。全磁気／affine 係数と全逆を保持する。

## 1. 戻すべき全 source

`TWO_QUADRICS_PURE_NIL_1_4_10_WHOLE_REDUCTION.md` §2 の literal 六括弧を使う。
元全十三係数と scalar／最高部／平方／線形因子の全補集合は旧独立 source 還元への依存である。
\(N=y\partial_x+z\partial_y+w\partial_z\) に固定し、全 \(N^2\)-normalizer 後の対象は
\(A=n_0\ne0,B=n_3,C=n_4\) 任意。
actual \(x'=x-2Cw\)（\(N^3\)-centralizer）で C を消すと、全六括弧は
\[
 P(A,B)=\left(
 Ah^2+Bhz-y^2/6,\ Bhw-yz/3,\ -yw/3-z^2/12,
 -z^2/4,\ -wz/2,\ -w^2/2\right).
\]
順序は \((xy,xz,xw,yz,yw,zw)\)、h 行は全零。
actual affine Pfaffian は \((2wy-z^2)^2/48-Ah^2w^2/2\)。
A 非零が相異二既約 quadric の正確な条件であり、B を零に限定していない。

## 2. 真の completed T20 と h-Casimir gauge

source は `reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T20.json` の全補正込みとする。
source 座標を \((a,b,c,d,H)=(x_0,x_1,x_2,x_3,x_4)\)、
period \((p_0,p_1)=(r,s)\)、smoothing \((t_0,t_1)=(t,u)\) とする。
真の raw source の \(-2tu\,bd/s\) correction を省かない。
source Euler 場を \(E_S\)、一次場を
\[
 L=(0,-rb,-(3r-s)c/2,-(r+s)d/2,0)
\]
とすると \(P_T-E_S\wedge L\) は H-Casimir であり、その六括弧は正確に
\[
\begin{array}{lll}
 P_{ab}=rab+td^2,&
 P_{ac}=(3r-s)ac/2-2tu\,bd/s,&
 P_{ad}=(r+s)ad/2,\\
 P_{bc}=(r+s)bc/2,&P_{bd}=-(r+s)bd/2,&P_{cd}=-rcd+ub^2.
\end{array}
\]
これは raw H 行の全四式から決めた Euler gauge であり、divisor だけを同定したものではない。
period open は \(rs\ne0\)。全 source の閉包が named T20 であることは
`research/notes/MPS_S61_CATALOGUE_20261103_v62.md` の既受理入力を使う。

実固定 frame
\[
 (h,x,y,z,w)=(H,a/6,d/2,b,c)
\]
の determinant は \(1/12\)、逆は \((a,b,c,d,H)=(6x,z,w,2y,h)\)。
以下 \(r=e,s=e+\sigma e^2,t=s,u=-s\) と置く。
generic \(e(1+\sigma e)\ne0\) で period open に入る。
\(-2tu/s=2s\) なので、この correction の period pole は binding 時点で正確に相殺する。

固定 frame 後の H-Casimir source は厳密に \(eK+\sigma e^2S\)。六括弧表は

| pair | K | S |
|---|---|---|
| xy | \(xy\) | \(xy/2\) |
| xz | \(xz+2y^2/3\) | \(2y^2/3\) |
| xw | \(xw+2yz/3\) | \(-xw/2+2yz/3\) |
| yz | \(yz\) | \(yz/2\) |
| yw | \(yw+z^2/2\) | \(z^2/2\) |
| zw | \(zw\) | \(zw/2\) |

## 3. 正確な nil flow と pushforward の符号

\(N=(y,z,w,0)\) の実 flow は
\[
 \phi_a(x,y,z,w)=
 (x+ay+a^2z/2+a^3w/6,\ y+az+a^2w/2,\ z+aw,\ w).
\]
determinant 一、逆は \(a\mapsto-a\)。actual pushforward の微分は
\(-\mathcal L_N\) であり、\(+\mathcal L_N\) ではない。
一次場の行列を \(D N\) とすれば \(\mathcal L_NP=N(P)-(DN)P-P(DN)^T\)。

\(P_0=P(0,0)\) とし、次の divergence零 bivector を置く。
\[
 J=(xz-y^2/3,\ xw+yz/3,\ (yw+z^2)/3,\ yw,\ wz,\ w^2).
\]
六式を微分するだけで
\[
 \mathcal L_NK=-2P_0,\quad\mathcal L_NS=J/2,
 \quad\mathcal L_NP_0=\mathcal L_NJ=0
\]
を得る。従って級数の打切りを仮定せず、**厳密に**
\[
 (\phi_{1/(2e)})_*(eK+\sigma e^2S)
       =P_0+eK+\sigma e^2S-\sigma eJ/4.                 \tag{1}
\]
右辺は全 e で多項式。これは有限 jet でなく実 tensor の全係数恒等式である。

## 4. 二つの actual affine translation と全六多項式

e 非零でさらに \(x'=x+dh/e\)、次に \(y'=y+3dh\) とする。
両方 determinant 一の実 PGL で、逆は同じ係数を引くこと。
\(P_0\) は x に依存せず、K／S／J は x について高々一次なので、
第一 translation の全負 e 冪は正確に消え、最終 tensor は全 e で多項式になる。

\(Y=y-3dh\)、\(A=3d^2/2\)、\(B=\sigma d/4\) と書けば、最終六式 \(Q_e\) は次の通り。
\[
\begin{aligned}
 Q_{xy}={}&Ah^2+Bhz-y^2/6
 +e(xY-\sigma dhY/2-\sigma xz/4+\sigma Y^2/12)+e^2\sigma xY/2,\\
 Q_{xz}={}&Bhw-yz/3
 +e(xz+2Y^2/3-\sigma xw/4-\sigma Yz/12)+2e^2\sigma Y^2/3,\\
 Q_{xw}={}&-yw/3-z^2/12
 +e(xw+2Yz/3+\sigma dhw/2-\sigma(Yw+z^2)/12)
 +e^2\sigma(-xw/2+2Yz/3),\\
 Q_{yz}={}&-z^2/4+e(Yz-\sigma Yw/4)+e^2\sigma Yz/2,\\
 Q_{yw}={}&-wz/2+e(Yw+z^2/2-\sigma wz/4)+e^2\sigma z^2/2,\\
 Q_{zw}={}&-w^2/2+e(zw-\sigma w^2/4)+e^2\sigma zw/2.
\end{aligned}                                                   \tag{2}
\]
h 行は全零。全 spatial degree は二、e degree も二以下。
generic e では、真 completed T20 source の h-Casimir lift に対する
固定 frame、\(\phi_{1/(2e)}\)、二 translation の actual pushforward そのものである。
逆順に三つの実逆を掛ければ元 source へ戻る。
全 frame determinant は固定 frame の \(1/12\) のみで、零にはならない。
raw T20 との違いは第2節の Euler 項を同じ actual frame で移したものだけ。

generic e の全 Jacobi は真 source から零、(2) が全多項式なので全 e で零である。
源の actual divisor は元全 tensor と同じ frame の signed Pfaffian 変換則で戻る。
特に e 零では厳密に
\[
                    Q_0=P(3d^2/2,\ \sigma d/4).
\]
その affine Pfaffian は \((2wy-z^2)^2/48-3d^2h^2w^2/4\)。
d 非零なら非零で、中心の projective tensor と対象 factor 型が保たれる。

## 5. Whole A非零・全B/C と元 source への結論

任意 \(A\ne0\)、任意 B に対し、\(d^2=2A/3\)、\(\sigma=4B/d\) を取る。
複素平方根の有限二選択で全 A を覆い、B 非零を反転する必要はない。
\(e(1+\sigma e)\ne0\) は中心を除く非空稠密開で、その全像は named T20。
(2) は中心でも非零の実多項式族なので、閉性から \(P(A,B)\) も T20 に入る。

C は第1節の \(N^3\)-変換の逆、元の \(n_5\) 等は既実 \(N^2\)-centralizer の逆で戻せる。
PGL5 の連結作用と非零 scalar は T20 を保つので、**型一 HQQ の全 actual source は T20**。
既 \(B=0\) の T39 境界認識とも矛盾せず、その部分が複数既知閉成分の交点に入る可能性を排除しない。

A 零は quadric 平方として本対象外、B 零も上の同じ曲線で保持した。
他の nil 型、型七の四未認識軌道、型十の \(\kappa\)-unit 軌道、全分類は別に残す。
本稿の新部分は literal whole curve と全実正逆である。
後続の有限 Laurent readback は同じ (1)–(2) と真 source だけを検証し、
別 anchor／operator／profile 又は higher jet を追加しない。
