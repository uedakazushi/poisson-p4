# 滑らかな quadric の可換 conormal は実 hyperplane 因子を強制する

2026-11-07。主義務 G、producer `plane_collision_review`。
**命題。** \(S\subset\mathbb P^4\) を smooth quadric surface とし、\(H=(h=0)\) をその唯一のspan hyperplaneとする。
\(\pi\) がSをPoisson invariantとし、\(I_S/I_S^2\) のconormal bracketが零なら、HもPoisson invariantである。
特に rank 四πの実退化quinticはhを因子に持つ。

**帰結。** 実既約被約quinticを持つrank四Poisson tensorの非正規曲面成分はsmooth quadricにならない。
restriction rank零・二の両方、零curve Rの全可約・非被約型、全normal-cone rankとvertical型を含む。
quadric cone又は他のSはこの線形非退化性を持たず、明示補集合である。
本稿は新しいpureproof producerであり、独立監査前。新CAS・source標本は無い。

## 1. projective tangent条件とEuler lineの正しいtwist

同次座標を \((h,x_1,x_2,x_3,x_4)\) とし、
\(I_S=(h,q(x))\)、\(q=x^{\mathsf T}Qx\)、Qはinvertible symmetric4×4 matrixとする。
Sのsmoothnessはこのrank四条件と同値。
Euler表示によりπを任意のhomogeneous quadratic raw bivector Πで表す。
全projective sectionがこの表示を持ち、gaugeは \(E\wedge W\)、
\(E=h\partial_h+\sum x_i\partial_{x_i}\)、W任意のlinear vector fieldである。

SのPoisson invarianceから、πのnormal dh,dq contractionはS上零。
rawの \(\{h,X_i\}_\Pi\)、\(\{q,X_i\}_\Pi\) は従ってEuler line上にある。
degreeを正確に保つと、あるglobal section

\[
 \lambda\in H^0(S,O_S(H)),\qquad\mu\in H^0(S,O_S(2H))
\]

により全iについて \(\{h,X_i\}|_S=\lambda X_i\)、\(\{q,X_i\}|_S=\mu X_i\)。
ここで \(X_0=h,X_i=x_i\)。
前者はtwisted Euler列のkernel \(O_S(H)\subset O_S(2H)^5\)、
後者は \(O_S(2H)\subset O_S(3H)^5\) に対応する。
局所rational coefficientをglobal sectionと無断で同一視したものではない。

quadricのrestriction maps
\(H^0(\mathbb P^3,O(1))\to H^0(S,O(H))\)、
\(H^0(\mathbb P^3,O(2))\to H^0(S,O(2H))\) はsurjective。
前者は四つのSegre coordinates、後者は九つのbidegree(2,2) monomialから直接確認出来る。
従ってλをxだけのlinear form、μをxだけのquadratic form \(x^{\mathsf T}Bx\) としてliftする。

## 2. 二つのnormal contractionを同時に消す全global gauge

linear field Wを

\[
 W(h)=\lambda(x),\qquad W(x)=\tfrac12 Q^{-1}Bx
\]

と取る。qがhに依らずQ invertibleなので、\(W(q)=\mu\) と \(W(h)=\lambda\) は独立に同時可解。
\(\Pi'=\Pi+E\wedge W\) とすれば

\[
 \{h,X_i\}_{E\wedge W}=hW(X_i)-X_iW(h),\qquad
 \{q,X_i\}_{E\wedge W}=2qW(X_i)-X_iW(q).
\]

よって両raw contractionはS上零となる。
したがって全iに対し \(\{h,X_i\}_{\Pi'},\{q,X_i\}_{\Pi'}\in I_S\)。
W自身がSにtangentである必要は無い。Euler gaugeはprojectiveπを全く変えない。
また \(E\wedge(\Pi')^2=E\wedge\Pi^2\) なのでactual退化quinticも変わらない。

## 3. intrinsic conormalからstrict normal bracketへ

任意のprojective chart \(X_j\ne0\) でnormal functionsを
\(u=h/X_j\)、\(v=q/X_j^2\) とする。
Π'が定めるprojective bracketを正確に展開すると

\[
 X_j^3\{u,v\}
 =\{h,q\}_{\Pi'}-\frac{2q}{X_j}\{h,X_j\}_{\Pi'}
 -\frac h{X_j}\{X_j,q\}_{\Pi'}.
\]

§2により両correctionはI2に入り、conormal可換性により左辺もI2。
従ってraw \(\{h,q\}_{\Pi'}\) はsheaf \(I_S^2\) に入る。
\((h,q)\) はprime regular sequenceなので、そのpowersはprimary・saturated。
実際 \(\operatorname{gr}_I R=(R/I)[H,Q]\) の各pieceはdomain R/I上freeであり、
g∉Iの乗法が単射となることからprimary性が従う。
ゆえに同次polynomialのordinary \(I^2\) membershipとしてもよい。
normal-line twistによるdenominator correctionをここで明示的に処理した。

## 4. 全raw h-rowと非退化quadricのradical

Π'のh-rowはdegree二でS上零だから、全i=1,…,4について

\[
 \{h,x_i\}_{\Pi'}=h\ell_i+q a_i,
 \qquad\ell_i\text{ linear},\ a_i\in\mathbb C.
\]

従って
\(\{h,q\}_{\Pi'}=h\sum_i\ell_iq_i+q\sum_i a_iq_i\)。
§3のI2 membershipをh0へrestrictすると \(q\sum_i a_iq_i\in(q^2)\)。
左辺はdegree三、q2はdegree四なので左辺は零。
polynomial ringはdomainなので \(\sum_i a_iq_i=0\)。
\(\nabla q=2Qx\)、Q invertibleから全 \(a_i=0\) である。

従ってΠ'の全h-rowはhで割れ、Hのprojective idealはPoisson invariant。
actual Pfaffian \(E\wedge(\Pi')^2\) の各項もhで割れる：
Eのh成分を使う項は明示的にh、他はΠ'のh-rowを一つ使うからである。
Euler gauge不変性から、これは元πのactual quinticのfactorでもある。
rank四・実既約被約quinticには矛盾する。

## 5. 非正規smooth quadric支持への全rankの適用

非正規支持のsurface成分SはHamiltonian flowで保存され、Poisson invariantである。
π|Sの一般rankは零又は二。
rank二ならS⊂SingDからactualPfaffianはI2に入り、接方向πがunitなdense openで
normal一次係数を比較して \(\{I,I\}\subset I^2\)。
smooth S上のbundle morphismなのでconormal bracketは全Sで零に延びる。
これは既受理 `IRREDUCIBLE_NONNORMAL_RANK2_POISSON_MODULE.md` のrank二補題そのもの。
rank零なら、既受理 `REDUCED_DEGENERACY_ZERO_SURFACE_ABELIAN_CONORMAL.md` と別作者reviewが、
被約actualDだけから全重複度conormal可換性を与える。

両rankとも本命題を適用してactualHfactorを得る。
Rのintegrality/reducedness、normal-cone root分裂、vertical factorの高次延長を追加仮定していない。
旧 integral-R exclusion と reducible-R root templates は必要条件として不変だが、
actual irreducible-D枝はこのglobal linear Euler argumentですべて排除される。
quadric coneではQにradicalがあり、最後のa_i零の段階が成立しないので別残余として保持する。
新しい成分数・他の支持・全G/Eの完了は主張しない。旧packet、global/main/Gitは変更しない。
