# canonical index二・degree五・rank二のcusp枝は空

2026-11-07、PRIMARY G。別作者による純証明。
対象は既約被約quintic Dのnormalizationがcanonical Gorenstein、`-KX=2L,L³=5`で、
導手の像Sがintegral非退化degree五、π|Sが一般rank二である枝だけである。
terminal、Q-factorial、smooth S、normal Γ又はreduced vertex curveを仮定しない。
新CAS又は数値検査は行っていない。凍結済みの長い≤24証明は変更しない。

## 全cusp sourceへの入口

既初回還元と新whole-source証明の§1–8から、cuspの導手divisorは`2Γ`、`Γ∼L`であり、
全closed inverse-residue formは

\[
\alpha=d(s_1/s_0),\qquad s_0,s_1\in H^0(X,L),\quad(s_0=0)=\Gamma
\]

と書ける。rank二からs1はs0に比例しない。導手quadratic qのpullbackは非零constant倍のs0²である。
導手schemeはCM pure、一般δ=1でreducedなので、integral支持Sそのものである。
qはrank三のirreducible quadraticである。

`W=ν*H0(P4,O1)`とする。s0がWに入ればΓの像Sがhyperplaneに入り、非退化性に反する。
s0,s1がmod Wで従属の場合は`s1=c*s0+ell`、ell∈Wとなる。
ell=0ならα=0でrank二に反する。ell≠0なら`q/ell²=(s0/ell)²`はπ|Dの非定数Casimirである。
`ell*dq-2q*dell`はEuler contraction零のΩ1(3)sectionで、そのHamiltonianはD上零。
従ってH0(TP4(-2))=0からambient上も零となり、一般rank四に反する。
q rank三なのでq/ell²は定数ではない。
これで**全従属marks**が排除され、独立marksをmissing generatorsに取れる場合だけが残る。
この排除でsource係数の稠密開だけを扱っているのではない。

## 全独立marksとliteralな二次恒等式

Y=s0、Z=s1をmissing degree一generatorに取る。
Beauville Theorem Bの証明は選択済みF0→Fのsurjectionを保持するので、この選択を保った対称表示を使える。
nonzero scalarを合わせ、全target座標のままで

\[
M=\begin{pmatrix}C&A&B\\A&u_0&u_1\\B&u_1&u_2\end{pmatrix},
\quad q=u_0u_2-u_1^2=Y^2,\quad \det M=D
\]

とする。u0,u1,u2は独立なlinear forms、A,Bは全quadratics、Cは全cubicsから始まる。
normalization algebraのadjugate pairingから

\[
Cu_2-B^2=qY^2=q^2.
\]

両辺はtargetのdegree四polynomialsであり、degree五のintegral D上の等式はliteral polynomial等式である。
u2で剰余を取ると`Bbar²=-u1⁴`。polynomial UFDで全場合は

\[
B=\varepsilon i u_1^2+u_2L,\quad\varepsilon=\pm1,\quad
C=(B^2+q^2)/u_2
\]

で尽くされる。Lは全五係数linear formであり、零係数も除外しない。
除算はliteral polynomialの割切れで、targetのu2≠0 chartだけへ制限していない。

\[
A_* = A-\varepsilon i u_0u_1-u_1L
\]

と置く。全式の展開から

\[
D=u_0q^2-2\varepsilon i u_1qA_*-u_2A_*^2,       \tag{1}
\]
\[
g_1=\operatorname{adj}(M)_{01}=qY
   =-\varepsilon i u_1q-u_2A_* .                 \tag{2}
\]

導手idealはmodule generators1,Y,Zへのqの乗法像なので

\[
I_S=(q,g_1,g_2),\qquad g_2=qZ=u_1A_*-qL.
\]

これは既whole-source表示の全cusp formulaであり、全extension coefficients、両sign、非terminal X、
非normal Γ、nonreduced vertex curveを保っている。

## degree五のintegral Sとの矛盾

qとg1はS上零である。Sは非退化なのでu2はその一般点で零ではない。
Sのintegralityと(2)から、AstarもS上恒等零となる。
ここではu2を全target又は全sourceで反転しておらず、integral支持上の一般点で積を消している。

rank三qはirreducibleである。Astarがqのscalar倍でなければ、二つのquadratics q,Astarは
regular sequenceとなり、そのproper complete intersectionはpure degree四である。
同じdimension二のintegral degree五surface Sを含むことは、degreeの加法性から不可能である。
従ってAstar=cq（c=0を含む）でなければならない。

しかし(1)はこのとき

\[
D=q^2\,(u_0-2\varepsilon i c u_1-c^2u_2)
\]

となる。u0,u1,u2の独立性から最後のlinear factorは非零であり、Dは被約既約でない。
仮定と矛盾する。

従属marksは第一節、独立marksの全sourceは第二・三節で処理したので、上記canonical index二の
**全cusp枝は空**である。特定のresidual rulingやgeneric cubic curveを前提にしない。
旧scroll/rank七/Gaの≤24上界は正しい保守的上界として残るが、この空性には不要である。
node枝及び他のnormalization・他のsupport degreeについては何も除外しない。
