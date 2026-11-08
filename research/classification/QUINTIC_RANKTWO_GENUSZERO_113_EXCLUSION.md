# 次数五 rank 二・種数零：全113型の排除と実射影の Schubert 補集合

2026-11-08、PRIMARY E、純幾何 producer proof。新しい算術は無い。

**命題。** `RANK_TWO_QUINTIC_SURFACE_GENUS_ZERO_SEPARATE_REVIEW.md` の
全113正規化 source は、actual 被約既約 quintic を持つ rank 四 Poisson tensor を与えない。
023／014 の全 source は本命題の明示的残余とし、同じ排除を流用しない。
また残る全113／023／014の実射影は \(\dim(V\cap K)=2\) を満たし、
旧上界に含めた intersection 三／四の Schubert strata は有限実射影では空である。

## 1. 受理済み全入口と正確な残余

旧 producer `RANK_TWO_QUINTIC_SURFACE_GENUS_ZERO_NONDOMINANCE.md` と
別作者 review の全源還元を使う。degree 五 integral 非正規支持 \(S\) は rank 二、
\(\delta_S=1\)、唯一の plane 残余 \(P\) も \(\delta_P=1\)。
全 generic cusp \(S\)、non-split node \(S\)、005、122、\(C_P=2F\)、014の \(2Q\) は既に排除されている。
ここで新たに扱う113では
\[
 Y=X=\mathbb P_{\mathbb P^1}(O(1)^2\oplus O(3)),\quad
 C_P=Q=\mathbb P(O(1)^2)\sim L-3F,
 \quad C_S=\Gamma_1+\Gamma_2\sim2L.
\]
\(\Gamma_i\) は異なる prime divisor、各々 \(|L|\) に属し \(S\) へ birational に写る。
元の finite birational map \(f:Y\to D\subset\mathbb P^4\) は任意の五次元 section subspace \(V\) による。
\(Q\simeq\mathbb P^1\times\mathbb P^1\)、\(L|_Q=O(1,1)\)、\(L^2Q=2\)。

旧 residue theorem は plane 上の residue を零とした。
generic node の逆形式・dualizing 形式の局所恒等式から、この零性は \(\pi|_P=0\) と同値である。
また global closed inverse-residue form は、\(\Gamma_i=(g_i)\)、\(g_i\in H^0(Y,L)\) として
\[
 \alpha=\lambda\,d\log(g_1/g_2),\qquad \lambda\ne0.               \tag{1}
\]
実際、この形は全 residues を持ち、差は全 divisor で正則な closed 一形式となり、
\(H^0(Y,\Omega_Y^1)=0\) から零。
これは選択した log ansatz のみでなく、旧全 closed form source の明示形である。

## 2. 実 node の rank 四条件は二つの正則 residue の非零 wedge

\(P\) の一般点で有限分離拡大・完成後、ambient smooth 座標を \(u,v,a,b\)、actual \(D=(uv=0)\) と取る。
actual Pfaffian は \(uv\) の unit 倍である。
\(\sigma=\pi^{-1}\) は closed、\(D\) に高々 simple pole を持つ。
SNC divisor に沿う closed simple-pole form は logarithmic なので
\[
 \sigma=A\,d\log u\wedge d\log v
  +d\log u\wedge(B\,da+C\,db)
  +d\log v\wedge(D\,da+E\,db)+F\,da\wedge db                    \tag{2}
\]
と書け、係数は正則である。
この log 性は \(d\sigma=0\) の二位極係数を見ることでも得られる：
\(u=0\) 上で normal \(du\) を含まない numerator は \(u\) で割れ、\(v\) についても同様。

逆行列の tangential 成分と Pfaffian の関係から、\(\pi|_P=0\) は \(A|_P=0\) を与える。
従って各 branch residue は \(P\) で正則で、その tangential restrictions は
\[
 \theta_+=(B\,da+C\,db)|_P,\qquad
 \theta_-=(D\,da+E\,db)|_P.
\]
他方 \(\sigma^2\) の \(d\log u\wedge d\log v\wedge da\wedge db\) 係数は
\(2(AF-BE+CD)\)。actual \(\operatorname{Pf}\pi=uv\times\mathrm{unit}\) なのでこの係数は unit であり、
\(A|_P=0\) より
\[
 \theta_+\wedge\theta_-\ne0.                                  \tag{3}
\]
normal units の選び方又は branch 交換はこの非零性を変えない。
係数体を Poisson 中心とは仮定していない。
non-split node でも有限分離 branch 拡大後の必要条件として (3) は有効である。

## 3. Q の全 double projection と dualizing residue の反不変性

\(f|_Q:Q\to P\) は有限 degree 二。
三つの \(O(1,1)\)-section による smooth quadric in \(\mathbb P^3\) の点からの線形射影で、
basepointfree なので中心は quadric 外にある。
全 off-quadric 点は quadric の automorphism で同じ型に移せる。
Segre matrix の中心を invertible antisymmetric matrix にすれば、その deck involution は
\[
 \iota(s,t)=(t,s)
\]
となる。従って \(\iota\) は二つの rulings を交換する。
任意の中心と section frame をこの同型の下で扱っており、一般中心を仮定していない。

quintic \(D\) の dualizing generator の rational pullback を \(\Omega_Y\) とする。
有限双対性により divisor は正確に \(-Q-\Gamma_1-\Gamma_2\)。
\(Q\) に沿う residue \(\beta=\operatorname{Res}_Q\Omega_Y\) は非零 rational 二形式で
\[
 \operatorname{div}(\beta)=-(\Gamma_1+\Gamma_2)|_Q.              \tag{4}
\]
これは adjunction と、\(\omega_Y(Q+\Gamma_1+\Gamma_2)\) の nowhere-zero trivialization の residue である。
marked curves の tangency／重複も右辺の scheme multiplicity のまま保持される。

generic node の同じ dualizing generator は、二 branch で
\(du/u\wedge\omega_P\) と \(-dv/v\wedge\omega_P\) の反対 residues を持つ。
従って function-field deck involution の下で
\[
 \iota^*\beta=-\beta.                                         \tag{5}
\]
ambient unit 自身が deck 固定とは仮定しない。同一の dualizing form の residue の functoriality を使う。
(4)–(5) から、\(a=g_1|_Q\)、\(b=g_2|_Q\) に対し
\[
 \operatorname{div}(ab)\quad\hbox{は }\iota\hbox{-不変}.         \tag{6}
\]
\(a,b\) は非零 \(O(1,1)\)-section。
もし一方が恒等的零なら prime \(\Gamma_i\) が \(Q\) を含み、異なる divisor class に反する。

## 4. 全 bilinear factor 補題

**補題。** \(a,b\in H^0(\mathbb P^1\times\mathbb P^1,O(1,1))\setminus\{0\}\) とし、
factor swap \(\iota\) が \(\operatorname{div}(ab)\) を保つとする。そのとき
\[
 \iota^*d\log(a/b)=d\log(a/b)\quad\hbox{又は}\quad-d\log(a/b).
                                                                    \tag{7}
\]
定数 ratio の場合も含む。

証明。bilinear polynomial は irreducible、又は二 rulings の linear factors の積である。
一つでも irreducible factor があれば、bihomogeneous polynomial 環の UFD により
\(\iota\) は \(a,b\) の irreducible \((1,1)\)-factors を固定又は交換する。
一方だけ irreducible ならそれを固定し、残りの二 linear factors の和も固定する。
従って ratio は非零定数倍の自身又は逆数となり、(7) が従う。

両方 reducible なら \(a=l(s)m(t)\)、\(b=n(s)o(t)\) と書ける。
(6) は binary quadratics \(ln\) と \(mo\) が比例することと同値。
\(l,n\) が独立なら、\((m,o)\) は定数倍を除き \((l,n)\) 又は \((n,l)\)。
前者では ratio は \(R(s)R(t)\)、後者では \(R(s)/R(t)\)、\(R=l/n\) であり (7) を与える。
\(l,n\) が比例すれば \(m,o\) も同じ factor に比例して ratio は定数となる。
共通 ruling、比例した marked sections、反復 factor を含め全場合を尽くした。証明終。

## 5. 全113型の矛盾

generic node の二 branch 上の \(\alpha\) の regular restrictions は、\(Q\) の function field で
\(\alpha_Q\) と \(\iota^*\alpha_Q\) に対応する。
(1) と (7) により
\[
 \alpha_Q\wedge\iota^*\alpha_Q=0.
\]
これは actual rank 四 node の必要条件 (3) に反する。
ゆえに全113 source は空。反復／reducible な marked conic、特殊射影、全 source frame を捨てていない。
有限個の点又は一般 mark の確認ではなく、全 UFD factor 補題を使った結論である。

## 6. 全三モデルの Schubert 補集合

旧全源では \(K=H^0(L-C_P)\) は四次元、\(\dim(V\cap K)\ge2\) とし、
intersection 二／三／四の上界13／9／3を保持した。
実際の finite image では \(C_P\) の像は全 plane \(P\) である。
従って restriction section 空間 \(V/(V\cap K)\) は少なくとも三次元。
他方二つの独立した ambient plane equations が \(V\cap K\) にあるので高々三次元。
ゆえに
\[
 \dim(V\cap K)=2.                                             \tag{8}
\]
これは113／023 の \(C_P=Q\)、014 の \(C_P=Q+F_b\) の全場合に有効。
023／014 の smooth scroll resolution \(Y\) の exceptional curve は small map で縮むが、
conductor surface の image dimension 二を変えない。
したがって intersection 三／四の補集合は、実 finite birational quintic image を与えず空である。
旧 dimension bound 自体を書き換えたり、その旧計算を false と呼んだりする必要は無い。

## 7. 残余と保存範囲

113 は全 actual rank 四で空。残る genus 零 rank 二 source は023／014の exact-intersection 二 stratum。
023 の quadric cone／その Hirzebruch resolution には factor-swap 補題を無条件で流用しない。
014 の \(Q+F_b\) は異なる二 prime を貼るので、同じ deck involution を仮定しない。
両者の全 \(\Gamma_i\)、plane restriction の零、singular／small 境界を残す。
全 E、既知親、新 component 又は他の genus／degree の結論は主張しない。
新 component の増分は零。新 CAS は無く、旧証明・global・主稿・Git は変更しない。
