# cusp-flex／cusp接線割線：全磁気を残したcubic因子上の非代数的葉

2026-11-08、PRIMARY E、`/root/log_manuscript`。全印刷sourceからの純手導出。
新しいpoint、operator、profile、minor、変形jet、CAS又は数値計算は無し。
全PGL parentの非所属、新成分性又は全分類の完成を本稿だけでは主張しない。

## 1. 加法対数留数による代数的leaf surfaceの障害

BをP1から有限点を除いた開曲線とし、tをその有理座標とする。
\(B\times\mathbb A^1_z\times\mathbb A^1_u\) 上の一般rank二foliationが

\[
 \eta=du+\omega(t),\qquad \omega\in\Omega^1_{\mathbb C(t)/\mathbb C}
\]

のkernelであるとする。P1上のomegaが一点でも非零留数を持てば、
Bの一般点を通る二次元の代数的leaf surfaceは存在しない。

実際、そのようなintegral surfaceはpartial_z方向に不変なので、
zを忘れた(t,u)-面内のintegral曲線Gammaの逆像となる。
これは局所UFDでも直接確認できる。surfaceのheight一primeをirreducible fで書くと、
partial_z fがfで割れ、z次数比較からpartial_z f零、従ってfはz非依存である。
一般leafではt方向が非零なのでGammaからBへの射は非定数。
projective closureのnormalizationを取ると、有限射 \(p:\widetilde\Gamma\to\mathbb P^1\) を得る。
uは同曲線上のmeromorphic functionで、foliationの接条件は

\[
 du=-p^*\omega.
\]

meromorphic functionのexact differentialの留数は全点で零である。
一方、pのramification indexをeとすればpullback留数は
\(e\operatorname{res}\omega\)、char零では非零のままであり矛盾する。
除いたpuncture上の点もcompact normalizationに戻るため、B上だけのalgebraic coverで逃れられない。
これはanalytic logarithmの形式表示だけでなく、任意の有限代数被覆を除く代数的証明である。

## 2. 第二cusp-flex源：全(I,J)非零で非代数的葉

全六磁気の正規化と実逆は
`TWO_LINES_PLANE_CUBIC_THREE_MAGNETIC_QUOTIENTS.md` §3及び独立
`TWO_LINES_PLANE_CUBIC_THREE_QUOTIENTS_SEPARATE_REVIEW.md` を使う。
ここではJを反転せず、全 \(R=Ih^2+Jhy\) を保持する。
全sourceはh-Casimir代表で

\[
 P=U\wedge\partial_z+V\wedge\partial_w+(-xw+Ih^2+Jhy)\partial_z\wedge\partial_w,
\]

\[
 U=2x(x+h)\partial_x+3xy\partial_y,\quad
 V=2xy\partial_x+3(y^2-(h+x)^2)\partial_y.
\]

actual cubic因子 \(C=h y^2-(x+h)^3\) のh一chartをnormalization parameter tで

\[
 x=t^2-1,\qquad y=t^3
\]

と書く。t非零の稠密開でこれはcubicのsmooth locusへのbirationalな実座標である。
非線形なこの正規化をPGL同値と呼ばない。
全restriction bracketは

\[
 \{t,z\}=A=t(t^2-1),\quad \{t,w\}=B=t^2(t^2-1)=tA,
 \quad \{z,w\}=K=-(t^2-1)w+I+Jt^3.
\]

t(t²-1)非零でrankは厳密に二。三変数(t,z,w)のkernel one-formは
\(Kdt-Bdz+Adw\)。A*tで割り

\[
 \eta=d(w/t-z)+\frac{I+Jt^3}{t^2(t^2-1)}dt
\]

を得る。\(u=w/t-z+I/t\) という全I,Jのrational座標（t単元上で逆も正則）を取れば

\[
 \boxed{\quad \eta=du+\frac{I+J}{2}\,d\log(t-1)
                       +\frac{J-I}{2}\,d\log(t+1).\quad} \tag{1}
\]

この式は部分fraction
\((I+Jt)/(t^2-1)=(I+J)/(2(t-1))+(J-I)/(2(t+1))\)
の全恒等式であり、source又は磁気parameterのgeneric値を代入した計算ではない。
(I,J)非零なら二留数が同時に零になることはない。
従って§1より、**その全点でcubic因子の一般点を通る代数的leaf surfaceは存在しない。**
特にJ単元のwhole lambda源(I,J)=(lambda,1)では全有限lambda、lambda正負一も含め同じ結論。
J零/I単元も直接(1)で覆う。I=J零は留数零であり、本障害から除外しない。
t=0,±1,無限遠のspecial fibersに別の代数的leafがある可能性はこの一般点の主張と区別する。

## 3. 第三cusp接線・尖点割線源：lambda非-1で非代数的葉

元全R6の正則shearとI単元chartのactual inverseは
`TWO_LINES_CUSP_TANGENT_SECANT_MAGNETIC_NORMAL_FORM.md` に固定される。
pair順(xy,xz,xw,yz,yw,zw)で全finite lambda源は

\[
 (0,-hx,3x(h-x),3(h-x)^2-3hy,-9xy,3xz+3\lambda h^2+3y^2).
\]

actual cubicは \(C_3=h^2y-(h-x)^3\)。h一でt=1-xとすればC3上でy=t³、
これは同chartの全cubicをA3(t,z,w)で表示する。
全restriction bracketは

\[
 \{t,z\}=1-t,\quad\{t,w\}=-3t(1-t),\quad
 \{z,w\}=3(1-t)z+3\lambda+3t^6.
\]

x=1-t非零でrankは二。kernelを(1-t)で割ると

\[
 \eta=d(w+3tz)+\frac{3(\lambda+t^6)}{1-t}dt.
\]

\[
 H_6(t)=t^6/2+3t^5/5+3t^4/4+t^3+3t^2/2+3t,
 \qquad u=w+3tz-H_6(t)
\]

と置く。H6'=3(t5+t4+t3+t2+t+1)だから

\[
 \boxed{\quad\eta=du-3(\lambda+1)d\log x,\qquad x=1-t.\quad} \tag{2}
\]

uへの変更はcubicのaffine chart上で全lambdaに正則なpolynomial座標であり、逆wも同じ次数の多項式。
projective又はPGLの変換と誤認していない。
lambda非-1ならx零で留数が非零。§1のfinite-cover residue論法より
**その全lambdaでcubic因子の一般点を通る代数的leaf surfaceは存在しない。**
lambda=-1ではeta=duとなり、本障害を適用しない。
元のI零磁気境界へgeneric陰性を外挿もしない。

## 4. Ambient leaf surfaceと、T37に使う際の正確な条件

cubic因子の一般点ではhと他のhyperplane因子が非零、cubic自体もsmoothであり、Poisson rankは二である。
したがってその点のambient symplectic leafはcubicに含まれる。
同点を通るintegralな二次元algebraic leaf surfaceがambient P4内にあれば、
そのrank二のleafは同surfaceの稠密開を張り、surfaceもcubic内に入る。
よって§2–3の障害はcubic内の特定の切片だけではなく、ambient algebraic leaf surfaceを除く。

**条件付き合成。** もしwhole T37に対して、actual Dの一般点を通るdegree有界な
algebraic invariant surfaceが存在するというproper closed Hilbert-incidence theoremを別に証明すれば、
上の非零留数源の全点はその条件を満たさずT37外となる。
generic T37 sourceだけのalgebraic integrabilityから、任意のsingular basis／degenerationへ無条件に外挿しない。
この全閉性・flat actualD・全componentの一般点を覆う量化が別途必要である。
本稿はそのparent側theoremの独立証明を先取りしない。

whole源閉包の親を除くには上記稠密開だけで十分だが、既知成分数の加算には
全prior成分及び他candidateとの比較と、全空間componenthoodの論理も別に必要。
それらを非代数的leafだけから宣言しない。
