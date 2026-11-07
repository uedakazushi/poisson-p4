# rank2 plane・smooth cubicのconormal接続補題：独立監査

2026-11-07。PLANE_SMOOTH_CUBIC_CONORMAL_CONNECTION.md の修正済み仮定を独立に確認し、
抽象flat-sl2補題とPoisson系への適用を受理する。
新しいsource計算・normal jet計算は行っていない。
receiptは CLASSIFICATION_PLANE_SMOOTH_CUBIC_CONORMAL_CONNECTION_REVIEW_20261107.json。

## 1. Poisson適用にはSが実Dの特異曲面である条件を使う

任意のrank2 invariant surfaceだけから \(\{I,I\}\subset I^2\) は従わない。
例えば局所
\(\pi=\partial_x\wedge\partial_y+u\partial_u\wedge\partial_v\)、
\(S=(u=v=0)\) では \(\{u,v\}=u\) である。
producer稿はこの指摘後、適用の仮定に **\(S\subset\operatorname{Sing}D\)** を明記した。
今回のPRIMARYである実非正規曲面成分ではこの条件を満たす。

必要なconormal bracket消滅は、Weinstein分解を使わず直接確認できる。
normal座標をu,v、接座標をx,yとする。
SのPoisson不変性によりmixed係数はIに入り、

\[
 \operatorname{Pf}\pi\equiv\pi_{xy}\pi_{uv}\pmod{I^2}.
\]

Sが滑らかな特異曲面なので \(D\in I^2\)。
\(\pi_{xy}\) がunitとなるdense openで \(\pi_{uv}\in I^2\) が従う。
bundle bracket \(\Lambda^2(I/I^2)\to I/I^2\) は正則なので、全Sへ零として延びる。

従って

\[
 \nabla_{df}s=\{\widetilde f,\widetilde s\}\bmod I^2
\]

はliftに依存せず、cotangent Lie algebroidのrepresentationとなる。
Jacobiがflatnessを与える。
このabelian conormal条件なしでは、lift変更によるinner bracketが残り、同じ推論はできない。

## 2. 一位極とhonestなglobal SL2 lift

S=P2上の \(\pi_S\) の零divisorをsmooth irreducible cubic Rとする。
U=S\Rではanchor \(\pi_S^\sharp\) が可逆。
局所では \(\pi_S=f\cdot(\text{非零正則二ベクトル})\)、R=(f=0)であり、
algebroid接続の正則係数へanchor inverseを掛けるので、ordinary connection matrixは
Rに沿う高々一位の極しか持たない。logarithmic性をここで仮定していない。

\(N^*=I/I^2=\mathcal O_{\mathbb P^2}(-1)\otimes\mathbb C^2\) である。
line factorのlocal frameを変えるtransitionはscalar \(g_{ij}I_2\)。
接続matrixの変更はscalarな \(d\log g_{ij}I_2\) だけなので、tracefree partは

\[
 \Omega\in H^0(\mathbb P^2,\Omega^1_{\mathbb P^2}(R))\otimes\mathfrak{sl}_2
\]

として全球にglueする。二つのnormal方向の固定basisはconstant GL2でしか変わらない。
局所matrixを \(A_i=\alpha_iI+\Omega\) と書くと、flatnessは
\(d\alpha_i I+d\Omega+\Omega\wedge\Omega=0\)。
traceから \(d\alpha_i=0\)、従って \(d\Omega+\Omega\wedge\Omega=0\) となる。

よって **d+OmegaはUの自明rank2 bundle上のgenuineなflat SL2接続** である。
projective monodromyだけをSL2へ未証明でliftしたのではない。
O(-1)のscalar flat holonomyが残る可能性は保持しており、元のN*接続そのものを零としない。

## 3. 非logならmeridianの固有値が−1となる

Rの局所方程式をx、接方向をyとして

\[
 \Omega=M(x,y)\frac{dx}{x}+N(x,y)\frac{dy}{x}
\]

と書く。M,Nはholomorphicでtracefree。
flatnessのx^-2係数は

\[
 [M_0,N_0]=N_0.
\]

nonlogなら一般点で \(N_0\ne0\)。
\([M_0,N_0^2]=2N_0^2\) のtraceから \(N_0^2=0\) であり、
\(N_0=E_{12}\) へ共役するとM0の二対角の差が1、下三角成分が零となる。
tracefree性からresidueの固有値は \(1/2,-1/2\)。

共鳴について追加仮定は不要である。
y固定のtransverse circle \(x=re^{i\theta}\) 上の平行移動は

\[
 \frac{dY}{d\theta}=-iM(re^{i\theta},y_0)Y.
\]

r→0で係数が一様にM0へ近づくので、matrix monodromyは
\(\exp(-2\pi iM_0)\) へ近づく。
異なるrでのmonodromyはradial平行移動によって共役なので、特性多項式はrに依存しない。
従って全rで両固有値が−1である。
これはresonant regular-singular系に対する通常の結論を直接導いており、
Jordan成分がある場合も固有値の議論は変わらない。

## 4. π1の外部入力をsmooth cubicへ正確に適用

Milnor, *Singular points of complex hypersurfaces* (1968), §6の設定と
Theorem6.5、printed pp.57–58を
[原著scanのtext](https://djvu.online/file/kCFrJq27uCPRi)で直接確認した。
孤立critical point、n>1という設定で、sphere fibrationのfiberがn-sphereのbouquetとなる定理である。
その証明全体を今回再証明したとはしない。

smooth projective cubic \(f_3=0\) では、homogeneous affine coneのcritical pointは原点だけ。
global fiber \(F=(f_3=1)\subset\mathbb C^3\) とsphere fiberの対応は

\[
 z\longmapsto\epsilon z/\|z\|,\qquad
 w\longmapsto f_3(w)^{-1/3}w
\]

である。後者ではsphere fiber上の正の実数 \(f_3(w)\) の正のcube rootを取る。
従ってn=2の定理が適用でき、Fはsimply connected。
F→P2\Rは自由なmu3作用のquotientなので \(\pi_1(U)=\mathbb Z/3\)。
局所meridianはcube-rootのdeck transformationのgeneratorへ写る。

独立照合として [Artal Bartolo–Dimca, arXiv:1507.08178v1](https://arxiv.org/pdf/1507.08178)
のIntroductionと§3、Proposition3.1の証明を直接読んだ。
そこではglobal Milnor fiberのcyclic covering、nodal曲線の場合の可換性、
および可換なcomplement群でのmeridian relationを区別して記している。
任意の特異既約cubicへ \(\pi_1=\mathbb Z/3\) を拡張していない。

genuine SL2 monodromyではmeridianの三乗がidentity。
固有値−1はこれに反するので、nonlogの場合は不可能である。

## 5. log形式消滅と結論の範囲

Omegaはlogarithmicとなる。
smooth irreducible Rのresidue完全列で、
\(H^0(\mathcal O_R)=\mathbb C\to H^1(\Omega^1_{\mathbb P^2})=\mathbb C\) は
1を \(c_1(R)=3H\ne0\) へ送る。
\(H^0(\Omega^1_{\mathbb P^2})=0\) と合わせて
\(H^0(\Omega^1_{\mathbb P^2}(\log R))=0\)、従ってOmega=0。

結論は、conormal接続が**projectively scalar**であること。
定数projective normal directionは全てhorizontalとなる。
二本の大域flat sectionがN*にある、元のscalar接続も零になる、
又はhigher normal項を含む全Poisson源を尽くした、という結論ではない。

singular/reducible/nonreducedなcubic R、SがSingDに入らない場合、
normal quadraticが平方又は零となる補集合は保持する。
blowup逆形式のvertical residueの高位極を、この別のconormal接続の一位極と混同しない。
このfirst-normal-jet補題だけからGの排除や全分類の完成は主張しない。
