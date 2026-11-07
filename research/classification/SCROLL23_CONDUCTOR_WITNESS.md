# S(2,3) の一つの有限射影：全導手の絶対既約性

2026-11-07。主義務 G、producer `plane_collision_review`。
これは **一つの処方済み射影の厳密証明書** であり、Poisson tensor の標本ではない。
相対導手族・留数・商の次元を使う G への矢印は
`RANK_TWO_SEXTIC_CONICAL_DIVISOR_NONDOMINANCE.md` に分けて証明する。
契約 `jobs/scroll23_conductor_witness.json` は gate PASS、実行時間は 0.358799 秒。
実行器は `code/classification_scroll23_conductor_witness.py`、不変の全データは
`reports/CLASSIFICATION_SCROLL23_CONDUCTOR_WITNESS_20261107.json` にある。
独立再生は未実施であり、以下は producer の証明と計算確認である。

## 1. 全域に定義された有限双有理射

\(Y=\mathbb P_{\mathbb P^1}(O(2)\oplus O(3))=\mathbb F_1\) とし、
tautological divisor を \(H\)、fiber を \(f\) とする。
\(H^2=5\)、\(Hf=1\)、\(K_Y=-2H+3f\)、負 section は \(C_0=H-3f\) である。
底座標 \((s:t)\)、fiber 座標 \((u:v)\) に対して四 section を

\[
 (Z_0,Z_1,Z_2,Z_3)=
 (us^2+vs^3,\;ust+vt^3,\;ut^2+vs^2t,\;v(s^3+t^3))
\]

とする。二つの係数行は
\(A=(s^2,st,t^2,0)\)、\(B=(s^3,t^3,s^2t,s^3+t^3)\) である。
全六 binary minor の gcd は 1 と厳密計算されたので、全 \(\mathbb P^1_\mathbb C\) で
二行は独立、従って四 section に共通零点はない。
この gcd 判定は有限体や実根だけの判定ではない。
手でも、\(s^3+t^3=0\) の場合だけ比例が問題になり、\(s=1\) で比例係数 1、
\(t=-1\) を強制した後、第三座標 \(t^2\ne t\) となる。
\(s=0\)、\(t=0\) はそれぞれ直接独立である。
四 section の線形独立性も全係数行列の rank 4 として確認した。

\(H\) は ample なので、これによる射 \(\nu:Y\to T\subset\mathbb P^3\) は finite。
四 section の独立性から像は plane でなく、
\(5=H^2=\deg(\nu)\deg(T)\) の素数性から \(\deg T=5\)、\(\deg\nu=1\) となる。
従って \(Y\) は実際の quintic surface \(T\) の正規化である。
以下で QQ の因数分解から選んだ五次式を、絶対既約と仮定してはいない。
像の integral 性と、この次数・双有理性が defining quintic の絶対既約性を与える。

## 2. 導手の canonical frame と全境界の復元

\(s=v=1\) で \(T=t/s\)、\(r=u/v\)、\(d=1+T^3\) と書くと

\[
 Z=(r+1,\ rT+T^3,\ rT^2+T,\ d),\qquad
 (x,y,z)=(Z_0,Z_1,Z_2)/d.
\]

\(d(y-Tx)+T-T^3\)、\(d(z-T^2x)-T+T^2\) の resultant から、
literal homogeneous pullback が零となる一意の五次因子 \(F(x,y,z)\) を取った。
全係数は report に保存され、式 \(F_{\rm hom}(Z)=0\) は厳密に確認された。
\(J=\det\partial(x,y)/\partial(T,r)\) に対し

\[
 d^3J=-T^3r-3T^3-3T^2-r=:J_{\rm num}.
\]

重要なのは、全導手の affine 表現が \(F_z/J\) だけではなく

\[
 G=d\,F_z(x(T,r),y(T,r),z(T,r))/J
   =F_{Z_2,{\rm hom}}(Z)/J_{\rm num}
 \tag{1}
\]

である点である。右辺の多項式除算の余りは零、非零 scalar は 1 だった。
得られた \(G\) は

\[
\begin{split}
G={}&-T^4r^2-3T^3r^2-2T^3r-3T^2r^3-3T^2r^2-2T^2r-T^2\\
 &-6Tr^3-6Tr^2+T-4r^3-6r^2-3r-1.
\end{split}
\]

式 (1) が全導手を与えることを、境界も含めて示す。
有限双対性より、smooth \(Y\) と Gorenstein hypersurface \(T\) に対し

\[
 \mathfrak c O_Y=\omega_Y\otimes\nu^*\omega_T^{-1},\qquad
 C=\nu^*K_T-K_Y=3H-3f.
\]

従って実導手は Cartier divisor であり、cycle の候補だけを計算しているのではない。
\(dx\wedge dy/F_z\) は \(\omega_T=O_T(1)\) の section \(Z_3\) に対応し、
pullback は \((J/F_z)dT\wedge dr\) となる。
\(F_\infty=(s=0)\)、\(F_d=(s^3+t^3=0)\) と書けば

\[
 \operatorname{div}(dT\wedge dr)=-2C_0-3F_\infty,
 \quad\operatorname{div}(d)=F_d-3F_\infty,
 \quad\nu^*(Z_3=0)=C_0+F_d.
\]

第一式は無限遠で \(T'=1/T\)、\(r'=r/T\) と置き
\(dT'\wedge dr'=-T^{-3}dT\wedge dr\) とすることで確認できる。
従って導手は

\[
 \operatorname{div}(F_z/J)+F_d+3C_0+3F_\infty
 =\operatorname{div}(G)+3C_0+6F_\infty.
\]

これは \(G\) の各項 \(r^iT^j\) を
\(u^iv^{3-i}s^{6-i-j}t^j\) とする \(3H-3f\) の正確な homogenization である。
全項で \(j\le6-i\) を確認した。得られた global section の二つの境界値は

\[
 G_{\rm glob}|_{s=0}=-t^4u^2v\ne0,\qquad
 G_{\rm glob}|_{v=0}=-su^3(4s^2+6st+3t^2)\ne0.
\]

従って chart \(sv\ne0\) の補集合に導手成分を落としていない。
\(d=0\) も式 (1) の polynomial identity により復元され、除外 fiber として残していない。

## 3. QQ の行列から複素数体上の絶対既約性へ

\(G\) は bidegree \((\deg_r,\deg_T)=(3,4)\) であり、
\(\gcd(G,G_r,G_T)=1\) を確認した。
未知多項式 \(A,B\) に次数制限

\[
 \deg_rA\le2,\quad\deg_TA\le4,\qquad
 \deg_rB\le3,\quad\deg_TB\le3
\]

を課し、閉一形式 \((A/G)dr+(B/G)dT\) の条件

\[
 G(A_T-B_r)-AG_T+BG_r=0 \tag{2}
\]

を全係数で線形化した。未知数は 31、非零係数行は 44。
全 \(44\times31\) 行列、全 row/column labels、核、minor index を report に保存した。
30 次 minor は \(3486368987873280\ne0\)、
\((A,B)=(G_r,G_T)\ne0\) が核にあるので、QQ でも C でも核の次元は正確に 1 である。

もし \(G=\prod_{i=1}^kG_i\) が C 上で相異なる絶対因子を持てば、各因子から

\[
 (A_i,B_i)=((G/G_i)(G_i)_r,\ (G/G_i)(G_i)_T)
\]

が上記の全次数制限と (2) を満たす。
対応する \(d\log G_i\) は \(G_j\) の一般点での residue が \(\delta_{ij}\) だから線形独立である。
従って \(k\le1\)。QQ での因数分解結果を C へ取り違えることなく、絶対既約性が従う。
§2 の境界チェックと合わせて、**この射影の全 global conductor は integral** である。

## 4. 証明の範囲

実行は処方済みの一射影だけで完了し、projection search、GB、有限体の因数分解、
Poisson 点の抽出を使用していない。
この計算から単独で全 Poisson の G 非優越や既知成分帰属は結論しない。
全射影族での導手の algebraic variation、split locus の properness、
固定退化因子上の留数の一意性、Aut の自由 fiber による次元の減算を次のノートで補う。
