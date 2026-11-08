# canonical index二 node 全源は T37：一つの正確な有理 PGL 同定

2026-11-08、PRIMARY E。全 source の純還元、一回の有限 QQ lift と literal な全括弧同定。
独立レビュー前の producer proof。新 component 数は零。

**定理。** 実退化因子が被約既約 quintic、integral 非退化 degree 五導手支持上の restriction が一般 rank 二で、
正規化が canonical Gorenstein index 二となる全 node source は既知成分 T37 に入る。
canonical／非terminal／nonreduced vertex curve、全 extension coefficients を保持する。
既済の全 cusp／従属 marks／single-prime 排除と合わせ、この canonical index 二枝の全 E 帰属を閉じる。

## 1. 全 source の補集合を戻す

`CANONICAL_INDEX_TWO_NODE_WHOLE_ORBIT_ENTRANCE.md` の全証明を使う。
旧全十一 parameter・全 \(i=0,1,2\) source のうち、\(i=1,2\) は conductor に plane を含むため空。
\(i=0\) でも二 extension forms の vertex quotient rank が零／一なら、
actual quintic と inverse-residue が未使用座標に依らず、horizontal affine brackets の次数三対五の比較から rank 四では空。
残る rank 二 extension は実 target coordinate basis を与え、非零 \(\kappa\) の fifth-root presentation change と完全な逆を持つ。
従って全 actual source の pair は
\[
 q=ac-b^2,\quad A=c^2+bd,\quad B=a^2+be,
 \quad C=c^2e+a^2d+2acb+bde-b^3,
\]
\[
 F=\det\begin{pmatrix}C&A&B\\A&a&b\\B&b&c\end{pmatrix},
 \quad g_1=bB-cA,\quad g_2=bA-aB,\quad\alpha=\lambda d\log(g_1/g_2),\ \lambda\ne0. \tag{1}
\]
へ移る。この還元は選択点への特殊化ではない。

## 2. 一回の全 ambient lift

gate 済み120秒／768MB契約 `canonical_index_two_node_canonical_ambient_lift.json` の一回だけを実行した。
runner は `code/classification_canonical_index_two_node.py`、stdlib Fraction のみ。
\(e=1\) chart における全六つの二形式恒等式
\[
 g_1g_2,i_\Pi\mathrm{vol}=dF\wedge(g_2dg_1-g_1dg_2)\pmod F
\]
を、一本の polynomial による通常除算と QQ 線形消去で解いた。
全150 raw coefficients、Euler24の冗長性と正しい divergence-zero gauge を保持し、全945行を保存した。
新しい Jacobi operator／GB／profile は作っていない。

実測0.028720秒、peak RSS13,696KiB。唯一の divergence-zero lift は次の \(P\) である。

|括弧|\(P_{ij}\)|
|---|---|
|ab|\(-6ab/5\)|
|ac|\(-2ac/5-2b^2\)|
|ad|\(2ad/5+4bc\)|
|ae|\(6ae/5-2bd-2c^2\)|
|bc|\(-6bc/5\)|
|bd|\(3bd/5+3c^2\)|
|be|\(-3be/5-3a^2\)|
|cd|\(-6cd/5+2a^2+2be\)|
|ce|\(-2ce/5-4ab\)|
|de|\(-de/5+2ac+3b^2\)|

全元線形行、divergence、元の350 projective Jacobi coefficients、126 signed actual Pfaffian coefficients を literal に検証し、
\[
 \operatorname{Pf}_{\rm actual}(P)=6F.                           \tag{2}
\]
全行順／列順／係数行列と producer bytes は report 内に保存した。
matrix SHA256 は `2baa36b2fb0da62fbdec44cec0802fa1a38b4c24af7395618c2d4b62a7e1892b`。

\(F\) の absolute irreducibility も弱い独立証明書で確認した。
\((d,e)\) の quadratic coefficients の gcd は \(b\)、constant term の \(b=0\) は \(-a^5-c^5\) なので content は一。
\((a,b,c)=(1,1,2)\) でその conic は
\(-2d^2+3de-e^2-13d+10e-22\)、homogenized matrix の determinant は \(1/4\)。
従って generic conic は代数閉包 \(\overline{\mathbb C(a,b,c)}\) 上も非退化、Gauss lemma から \(F\) は absolute irreducible。
この空間点での非零 determinant は polynomial 非零性の証明書であり、Poisson source の標本ではない。
\(g_i\) は非零 cubic なので \(D\) 上で恒等零でなく、dense spatial chart の全補集合は復元される。

## 3. 真の completed T37 と有理行列

source は唯一、`reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T37.json` の
\(p_0=t_0=t_1=t_2=t_3=1\) とする。未完成の naive T37 は使わない。
source coordinates を \(X=(x_0,x_1,x_2,x_3,x_4)\)、target を \(Y=(a,b,c,d,e)\) として
\[
 Y=M X=(2x_3,\ x_0,\ -2x_4,\ -12x_2,\ 12x_1).                  \tag{3}
\]
\(\det M=576\)、逆は
\((x_0,x_1,x_2,x_3,x_4)=(b,e/12,-d/12,a/2,-c/2)\)。
source から target への raw pushforward は、全十括弧を代入すると次の通り。

|括弧|\((M_*P_{T37})_{ij}\)|
|---|---|
|ab|0|
|ac|\(-2ac+2b^2\)|
|ad|\(-ad-4bc\)|
|ae|\(-3ae+2bd+2c^2\)|
|bc|0|
|bd|\(-3c^2\)|
|be|\(3a^2\)|
|cd|\(3cd-2a^2-2be\)|
|ce|\(ce+4ab\)|
|de|\(-de-2ac-3b^2\)|

これは canonical table と比べて正確に
\[
 M_*P_{T37}=-P+E_Y\wedge R,\qquad
 R=(6a,0,-6c,3d,-3e)/5.                                      \tag{4}
\]
全非零係数を表に記載し、他は零なので、表は全150 raw coefficients の恒等式である。
source の divergence は \((0,-3x_1,3x_2,6x_3,-6x_4)\) であり、
その canonical Euler gauge を取れば (4) は単に \(M_*P_{T37,\mathrm{div}0}=-P\) となる。
逆向きも (3) の実逆行列と (4) から直接従う。

signed actual Pfaffian の変換則と (2) は
\[
 576F_{T37}(X)=6F(MX),\qquad F(MX)=96F_{T37}(X).                 \tag{5}
\]
Euler 項は projective Pfaffian に寄与せず、tensor scalar \(-1\) の二乗は一。
ここで tensor／Pfaffian の source と target の向きを固定した。

## 4. 全 actual source への矢印

任意の元の canonical node actual tensor を第1節の可逆変換で (1) へ移す。
同じ \(D,\alpha\) に対応する projective ambient tensor は高々一つである。
これは dense smooth \(D\) 上の \(\alpha=i_{\pi_D}\Omega_D\) と
\(H^0(\Lambda^2T_{\mathbb P^4}(-5))=0\) による旧受理済みの全 source 一意性。
第2節の \(P\) がその tensor で、第3節により実際に T37 と同型である。
従って**全対象点が T37 の PGL 軌道に入る**。
残る source が canonical normalization を実際に与えることを必要条件の表だけから逆向きに仮定していない。
この包含にはその converse は不要である。

全 dependent marks／single-prime／non-split node と cusp は旧全証明で空、
独立 node の全三UFD型／全十一parameter／低 vertex rank は第1節に戻している。
したがって canonical index 二枝の全 E は T37 又は空で覆われる。
正の一標本から全支持を推論したのではなく、先に全 source の可逆還元と補集合を証明した。

新 component、T37 の全 converse、全 E、他 normalization／degree／factor の主張はしない。
一回の finite lift 以外の追加 ambient operator は無く、後半は実有理行列の全 literal identity である。
旧 reports／global／主稿／Git は変更せず、同名 manifest に exact inputs を固定する。
