# canonical index二 node 全源：一つの実軌道又は空への入口

2026-11-08、PRIMARY E、純幾何 producer partial theorem。独立レビュー前。
新しい行列計算、GB、minor、jet 又は source 点列は実行しない。

**到達点。** 既受理の canonical index 二 node の全十一 parameter sources を戻すと、
実 rank 四 source は一つの canonical \((D,d\log(g_1/g_2))\) pair に移る。
その pair に ambient tensor が存在すれば、全対象は一つの actual PGL 軌道である。
存在・rank 四・既知親はまだこのノートから主張しない。

## 1. 全 source と既に排除した補集合

`RANK_TWO_QUINTIC_CANONICAL_INDEX_TWO_NONDOMINANCE.md` §§1–9、
その別作者 review、及び全 cusp 排除を入力とする。
actual \(D\) は被約既約 quintic、\(S\) は integral 非退化 degree 五支持で一般 rank 二。
正規化 \(X\) は canonical Gorenstein、\(-K_X=2L\)、\(L^3=5\)。
terminal／Q-factorial 又は smooth marked hyperplane を仮定しない。

single-prime／non-split node、cusp、mod ambient sections で従属な marks は全て既に排除されている。
残る独立 node marks \(Y,Z\) と unit を保持した全表示は
\[
 M=\begin{pmatrix}C&A&B\\A&u_0&u_1\\B&u_1&u_2\end{pmatrix},
 \quad q=u_0u_2-u_1^2=YZ,
 \quad F=\det M,
 \quad AB-u_1C=q^2.                                             \tag{1}
\]
\(u_0,u_1,u_2\) は独立、\(q\) は rank 三。
全 conductor ideal は
\[
 I_S=(q,g_1,g_2),\quad
 g_1=u_1B-u_2A=qY,\qquad g_2=u_1A-u_0B=qZ.                     \tag{2}
\]
これは integral \(S\) 自身の ideal であり、generic support だけの同定ではない。
全 closed inverse-residue form は
\(\alpha=\lambda d\log(Y/Z)=\lambda d\log(g_1/g_2)\)、\(\lambda\ne0\)。

旧 UFD source は全 \(i=0,1,2\) に対し
\[
 A=\kappa u_0^i u_2^{2-i}+u_1L_1,
 \quad B=\kappa^{-1}u_0^{2-i}u_2^i+u_1L_2,
 \quad C=(AB-q^2)/u_1,\quad \kappa\in\mathbb C^*,               \tag{3}
\]
\(L_1,L_2\) は全 target linear forms である。
除算は literal polynomial divisibility であり、\(u_1\ne0\) の空間 chart を選んだのではない。
以下、この全零係数・全 vertex 方向を保持する。

## 2. i一／二は全係数で conductor に plane を含む

plane \(\Pi=(u_0=u_1=0)\subset\mathbb P^4\) を使う。
\(i=1\) なら \(A|_\Pi=B|_\Pi=0\)、\(i=2\) でも \(A|_\Pi=0\)。
従って両場合で (2) の全三 generator は \(\Pi\) 上零である。
これは \(L_1,L_2\) の全値で成立する。
しかし integral degree 五 surface \(S\) が同じ dimension 二の plane を含むなら、
その reduced support は plane 自身となって degree に反する。
よって全 \(i=1,2\) は空。
残る唯一の source は
\[
 A=\kappa u_2^2+u_1L_1,\qquad
 B=\kappa^{-1}u_0^2+u_1L_2.                                   \tag{4}
\]

## 3. basic inverse-residue を持つ quintic cone の補題

**補題。** \(\mathbb P^4\) の rank 四 Poisson tensor が actual 被約既約 quintic
\(F(x_0,x_1,x_2,x_3)=0\) を持ち、\(F\) は第五座標 \(w\) に依らないとする。
同じ標準 dualizing generator による inverse-residue 一形式 \(\alpha_D\) が
\(w\) を含まない base surface から来るなら矛盾する。

証明。\(x_0=1\) の chart を取り、base coordinates を \(x_1,x_2,x_3\) とする。
\(F\) は irreducible degree 五で linear factor を持たないので、その dehomogenization も degree 五。
smooth な base surface 上の局所座標を \(a,b\) と書けば
\(\Omega_D=dw\wedge\Omega_{\rm surface}\)（非零定数 convention を除く）。
\(i_{\pi_D}\Omega_D=\alpha_D\) に \(dw\) 成分が無いので、
\(\pi_D\) の base surface への horizontal bivector は零である。
従って三つの base coordinate pairs \(\{x_i,x_j\}\)、\(1\le i<j\le3\)、は \(D\) 上零。

任意の global projective quadratic bivector のこの chart 上の三 bracket は degree 三以下の polynomials。
具体的には homogeneous raw quadratic \(P\) に対し
\(\{x_i,x_j\}=P_{ij}-x_iP_{0j}+x_jP_{0i}\) である。
irreducible degree 五の \(F\) 上零の degree 三以下 polynomial は零多項式なので、
全三 base brackets は ambient でも零。
すると affine bivector は \(\partial_w\wedge V\) の形で rank 二以下となり、rank 四に反する。証明終。

これは cone そのものを無条件に禁止する補題ではない。
**同じ actual dualizing form に対する residue が vertex 座標を含まない**ことを必要とし、以下でそれを literal に確認する。

## 4. 全 vertex-rank 零／一の排除

\(U=\operatorname{span}(u_0,u_1,u_2)\subset H^0(O_{\mathbb P^4}(1))\) とし、
\(r=\dim\operatorname{span}(L_1,L_2)\pmod U\) を取る。
\(r\le1\) なら \(u_0,u_1,u_2,L_1,L_2\) は高々四つの独立 target coordinates しか使わない。
残りの一つを \(w\) と選ぶ。
(1)、(2)、(4) から \(F\)、\(g_1\)、\(g_2\) はいずれも \(w\) に依らない。
標準 quintic dualizing generator は affine residue volume なので、上の cone 補題と同じ \(dw\) factor を持つ。
\(\alpha=\lambda d\log(g_1/g_2)\) は literal に \(w\)-basic である。
従って補題が全 \(r=0,1\) を排除する。
vertex coefficients の一つを generic と仮定したのではなく、両方の低 rank 補集合を処理した。

残る全 actual source は \(r=2\)、すなわち
\[
 (u_0,u_1,u_2,L_1,L_2)\quad\hbox{が target の実 linear coordinate basis}      \tag{5}
\]
である。

## 5. 全逆を持つ一つの canonical pair

旧 source の実 presentation change を使う。
\[
 Y'=tY,\ Z'=t^{-1}Z,\quad
 u_0'=u_0/t^2,\ u_1'=u_1,\ u_2'=t^2u_2,
 \quad A'=A/t,\ B'=tB,\ C'=C.                                 \tag{6}
\]
これは symmetric congruence で \(\det M\)、\(q\)、unit/algebra と \(d\log(Y/Z)\) を保存する。
\(L_1'=L_1/t\)、\(L_2'=tL_2\)、\(\kappa'=\kappa/t^5\)。
\(t^5=\kappa\) を選べば \(\kappa'=1\)。これは有限五通りの選択で、\(\kappa=0\) は元の UFD source に無い。
逆は \(t\) を \(t^{-1}\) に替えた同じ式であり、target の空間点を反転していない。

(5) から target coordinate basis を
\((a,b,c,d,e)=(u_0',u_1',u_2',L_1',L_2')\) と出来る。
全 actual \((D,\alpha)\) はこれにより
\[
 q=ac-b^2,\quad A=c^2+bd,\quad B=a^2+be,
 \quad C=c^2e+a^2d+2acb+bde-b^3,                              \tag{7}
\]
\[
 F=\det\begin{pmatrix}C&A&B\\A&a&b\\B&b&c\end{pmatrix},\quad
 g_1=bB-cA,\quad g_2=bA-aB,\quad\alpha=\lambda d\log(g_1/g_2).\tag{8}
\]
\(C=(AB-q^2)/b\) は展開した polynomial (7) であり、\(b=0\) chart を失わない。
元の全係数と完全な target basis は (6) と (5) の逆で戻せる。
すなわちこれは選択点の例ではなく、全 actual source の可逆正規化である。

## 6. 未解決の入口を正確に一つにする

旧 restriction injectivity
\(H^0(\mathbb P^4,\Lambda^2T(-5))=H^0(\Omega^2)=0\) と、
\(\alpha=i_{\pi_D}\Omega_D\) の injectivity により、(8) に対応する projective ambient tensor は高々一つ。
従って全 canonical index 二 rank 二 source は、**この tensor が存在すれば一つの actual PGL 軌道、存在しなければ空**となる。
ここではまだ canonical pair の ambient lift／rank 四／normalization 全条件を算術で確認していない。
実際の答えが既知親の境界か新候補かも未判定である。

次の弱い有限課題は、(8) の一つの全 pair に対する ambient linear lift を exact に調べること。
例えば affine volume \(\mathrm{vol}\) に対し
\[
 g_1g_2, i_\Pi\mathrm{vol}
  =dF\wedge(g_2dg_1-g_1dg_2)\pmod F                           \tag{9}
\]
という全 restriction 方程式は \(\Pi\) の150 raw coefficients に線形である。
閉じた \(\alpha\) のため、もし正しい ambient lift が存在すれば、Jacobi は
\(H^0(\Lambda^3T(-5))=0\) から自動的に従う。
actual Pfaffian が非零なら定数倍の \(F\)、零ならこの canonical pair は rank 四 source を与えない。
新しい parameter profile、点 farm、GB 又は長い jet は必要ない。
この finite linear task は別の具体的 job 契約・review 後だけに実行する。

新 component 数、全 canonical branch の pointwise membership、全 E 又は総分類は未主張。
今回保存するのは全補集合を戻した純粋な入口であり、旧 whole-source proof／台帳／主稿／Git を変更しない。
