# 零symbol全39源の一実軌道：別作者による純証明監査

2026-11-08、producer root、reviewer plane_collision_review。
`PLANE_ZERO_SYMBOL_WHOLE_SINGLE_ORBIT.md`（SHA256 `150f8085d666314e2d65100d69c5e73d61e285b016c85c8f7d5c21d1167516b2`）を
全読解し、**全対象を一つの明示 actual PGL 軌道へ移す新しい二つの矢印を PASS とする**。
新 CAS／operator／rank／minor／kernel／点の計算は無い。
ambient 成分性・既知親・新 component 数は未解決のままである。

## 1. 再利用する旧全源と今回の独立範囲

旧 `PLANE_RANK_ZERO_ZERO_SYMBOL_LIE_ALGEBROID_BOUND.md` とその別作者 review は受理済み入力とする。
旧全39源・二chart・全weightの bridge を新たに再計算したとは言わない。
その全結論は、対象の \(\Pi=\Pi_1+\Pi_2\) について
\[
 A=O_{\mathbb P^1}(-1)^3\simeq\operatorname{Der}_{\rm Lie}(K),
 \quad0\to L=O(-5)\to K\xrightarrow\chi O\to0,
 \quad e\in H^1(L)\setminus\{0\},
\]
\(K\) 全fiber aff1、\(\Pi_2\) は全九係数の global CE 二cocycleというもの。
anchor共通零／恒等零、abelian \(K\)、split \(e=0\) の全補集合は旧定理で因子又はrank二以下へ排除されている。
今回独立に確認したのは repeated-root extension の全排除、全二cocycle exact性、
それらを actual linear maps へ戻す矢印と明示 canonical model である。

## 2. infinitesimal projective stabilizer の全 lift

trace \(\chi\) と derived line \(L\) は Lie bundle \(K\) 自身から canonical に回復される。
base identity 上の Lie automorphism は quotient \(O\) を固定し、\(L\) を constant scale \(c\) で移し、
残りは \(\operatorname{Hom}(O,L)=H^0(L)=0\)。
extension compatibility は \(ce=e\)、非零 \(e\) より \(c=1\)。
従って base identity 上の infinitesimal kernel は零である。

逆に \(\theta\in\mathfrak{sl}_2\) が \([e]\) の projective stabilizer に属すれば、
\(\theta e=\lambda e\) と書ける。
\(O(-5)\) は SL2 の usual linearization を持つ。
dual numbers 上でその base action と \(L\) の \(-\lambda\) infinitesimal scale を合成すると、
extension class 自身の variation は零。
Čech cocycle の差は coboundary となるので、identity を中心とする extension 同型が得られる。
これは quotient \(\chi\) を保ち、
\([k_1,k_2]=\chi(k_1)k_2-\chi(k_2)k_1\) も保つため、Lie derivation の global section を与える。
単に underlying vector bundle の lift を Lie lift と読み替えていない。
odd degree の line bundle を PGL2-linearize せず、SL2 と line scale を使ったことも重要である。

したがって
\[
 H^0(\operatorname{Der}_{\rm Lie}K)
   \simeq\operatorname{LieStab}_{\operatorname{PGL}_2}([e]).
\]
Serre duality の binary cubic 表現で、double＋simple root と triple root の stabilizer はそれぞれ正次元。
左辺は \(H^0(O(-1)^3)=0\) なので両型は空。
三相異根型は一つの PGL2 軌道で、全 source の \(e\) はこの型に属す。
全 underlying \(A\) frame の変更は constant GL3、base action は normal GL2 lift。
いずれも \(\mathbb P^4\) の actual block linear map であり、抽象同型だけに止まっていない。

## 3. transgression は Atiyah line bundle class に一致する

derived line \(L=[K,K]\) は全 Lie derivation に保存される。
restriction は自然な map
\[
 A=\operatorname{Der}_{\rm Lie}K\longrightarrow\operatorname{At}(L)
\]
を与える。symbol は同じ base vector field、inner \(\operatorname{ad}k\) の \(L\) 上の作用は
\(\chi(k)\operatorname{Id}_L\)。従って kernel は正確に \(L\)、rank 二 quotient は
\[
 A/L\simeq\operatorname{At}(L).
\]
局所的な symbol と scalar multiplication の全 lift があるので surjectivity も確認出来る。
よって trace character による transitive extension の pushout は line bundle の Atiyah extension そのもの。

旧 analytic CE の \(\mathcal H^1=\mathbb C\chi\) はこれと自然で、
Hochschild–Serre の
\(d_2:H^0(\mathbb P^1,\mathbb C\chi)\to H^2(\mathbb P^1,\mathbb C)\)
は \(c_1(L)\)（sign／\(2\pi i\) convention を除く）。
Čech horizontal lifts の差を \(\chi\) へ写すと \(d\log\) of transition of \(L\) になるためである。
\(\deg L=-5\) なのでこの map は非零、両辺一次元より同型。
holomorphic connection が存在すると置いて curvature を零にしていない。

局所 CE cohomology は q零／一だけ、base constant cohomology は p零／二だけなので、
この differential は total degree 一／二の全項を消す。
従って全 CE hypercohomology は \(H^1=H^2=0\)（参考として \(H^0=H^3=\mathbb C\)）。
各 term \(O,O(1)^3,O(2)^3,O(3)\) は analytic higher cohomology が零なので、
これは**同じ全 global cochain complex**の cohomology に等しい。
algebraic global sections と differential も同じであり、非 O-linear differential に対する GAGA を仮定していない。

全九係数の二cocycle は従って全て \(d_A\theta\)、\(\theta\in H^0(A^*)\simeq\mathbb C^6\)。
\(H^1=0\) と global \(d_0=0\) より \(d_A\) はこの六次元空間上 injective。
旧≤七という上界の一つの cohomology class を本当に消した新しい結論である。

## 4. 全 cocycles を actual finite shear で消す

fiber coordinates \(x_i\)、base \(t\) で \(y_i=x_i+\theta_i(t)\) とすると、
linear bracket の constant 部分は
\[
 c'_{ij}=c_{ij}+\rho_i\theta_j-\rho_j\theta_i-\sum_k f_{ij}^{k}\theta_k
        =c_{ij}+(d_A\theta)_{ij}.
\]
\(\{\theta_i(t),\theta_j(t)\}=0\)、元 bracket は fiber-linear なので高次の有限作用項は無い。
\(\theta_i\in H^0(O(1))\) だからこれは全域の実 linear shear
\(a_i\mapsto a_i+\alpha_iw+\beta_i z\)。逆は符号反転。
これで全 \(\Pi_2\) を消し、前節の唯一の \(\Pi_1\) actual block orbit へ戻る。
generic cocycle 又は formal shear だけを扱ったものではない。

## 5. explicit model の独立手確認

\(b(t)=t^{-1}+t^{-4}\)、\(g=t^{-5}\) の producer transition を微分して
\[
 d_\infty=-t^2d-5tH-(t^2b'+5tb)E,
 \quad H_\infty=H+bE,\quad E_\infty=gE
\]
を独立に得た。
producer の \(\operatorname{Der}(K)(1)\) 三sections は finite chart で determinant 三。
\(r=1/t\) で O(1) transition も含めた反対 chart の三列は
\[
 (-1-r^3,-5r^2,4r),\quad(-r^2,-4r,3),\quad(-r,-1,0),
\]
全て regular、determinant は非零定数 \(-3\)。
従って有限 chart だけの trivialization ではなく、全 bundle が確かに \(O(-1)^3\) となる。

局所 \([H,E]=E\) と symbol 微分から三つの表示 bracket を全成分で再導出した。
raw homogeneous tensor を w一 chart へ戻すと、anchor は
\((1+t^3,t,t^2)\)、fiber brackets は
\[
 \{a,b\}=a+t^2b-4tc,\quad
 \{a,c\}=-2ta+4b+t^2c,\quad\{b,c\}=c,
\]
であり、producer の全 ten raw brackets と一致する。
従ってその全 projective Jacobi は Lie algebroid の Jacobi から従う。

signed Euler-Pfaffian を五つの four-by-four terms に分けて手展開すると、正確に
\[
 F=bz^4-4cz^3w+3az^2w^2-4bzw^3+cw^4
\]
を得た。係数 convention に依る「非零らしい」確認ではない。
base 三係数の binary quartics の gcd は一、primitive linear polynomial なので絶対既約。
normal 二方向の定数方向微分は a係数から零と強制され、残る三 quartics は独立なので essential 五変数。
base 三偏微分の共通零は \(z=w=0\) のみで、Sing(D) の reduced support はその唯一の plane。
zero anchor／abelian／split／repeated cubic の全補集合を処理した後の実 rank 四 model である。

## 6. 独立 verdict と未解決の親

全対象の \(e\) は一つの actual base orbit、全 \(A\) frame は actual GL3、
全 \(\Pi_2\) は actual finite six-shear で消える。
従って全39源の actual rank 四・被約既約 quintic 部分が一つの明示 actual PGL 軌道へ入るという命題を受理する。
これが新たに確認した E の構造還元である。

その軌道の ambient tangent／成分性／既知成分への包含又は新 count は本査読で確認していない。
新 CAS は無く、producer 全 bytes、旧 source／global／主稿／Git は変更しない。
exact 入力と review receipt は同名 manifest に保存する。
