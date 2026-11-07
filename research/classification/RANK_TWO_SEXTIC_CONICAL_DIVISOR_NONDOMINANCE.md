# rank 二・次数六支持：conical quintic の全像は次元二十四以下

2026-11-07。主義務 G、producer `plane_collision_review`。
**命題。** \(\pi\) を \(\mathbb P^4_\mathbb C\) 上の rank 四 Poisson tensor とし、
実退化因子 \(D\) を既約被約 quintic とする。
非正規支持に integral surface \(S\) があり、\(\deg S=6\)、\(\pi|_S\) の一般 rank が二、
かつ \(D\) の正規化が singular minimal-degree scroll である場合、
この全 projective Poisson 像の次元は高々 24 である。
新成分・既知成分帰属・点の全被覆は主張しない。

本ノートは producer の幾何証明であり、独立監査は未実施である。
唯一の新計算は `SCROLL23_CONDUCTOR_WITNESS.md` の一射影で、
その全厳密証明書は `reports/CLASSIFICATION_SCROLL23_CONDUCTOR_WITNESS_20261107.json`。
smooth 二型を扱う別作者ノート
`RANK_TWO_SEXTIC_SURFACE_NONCONICAL_DIVISOR_EXCLUSION.md` の結論は、本命題の入力にはしない。

## 1. rank に依存しない既受理の正規化入口

`RANK_ZERO_SEXTIC_CONDUCTOR_EXCLUSION.md` §§1–2 とその別作者監査
`RANK_ZERO_SEXTIC_CONDUCTOR_SEPARATE_REVIEW.md` §§2–3 の論証は、restriction rank を使わない。
従って現在の仮定でも、\(\nu:X\to D\)、\(L=\nu^*O_D(1)\) に対して

\[
 \delta_S=1,\quad g(X,L)=0,\quad L^3=5,\quad h^0(X,L)=8
\]

となり、二次元の非正規支持は \(S\) 一つだけである。
完全 \(|L|\) は \(X\) を \(\mathbb P^7\) の degree 五 minimal-degree threefold として埋め込み、
\(\nu\) はその全 basepoint-free 五 section の線形射影である。
全五型のうち、本ノートで扱う singular 三型は

\[
 X_{014}=S(0,1,4),\qquad X_{023}=S(0,2,3),\qquad X_{005}=S(0,0,5).
\]

rank 零専用の cone 排除は rank 二の本件へ流用しない。
最初の二型の smooth scroll resolution
\(p:\widetilde X=\mathbb P(O\oplus O(a)\oplus O(b))\to X\) は
zero-summand section \(\mathbb P^1\) のみを vertex 点へ縮める。
従ってこれは small resolution で、\(X_{\rm reg}\) 上同型、exceptional locus は codimension 二。
\(\operatorname{Cl}X\simeq\operatorname{Pic}\widetilde X=\mathbb ZL\oplus\mathbb ZF\)、
\(K_{\widetilde X}=-3L+3F\)、導手の strict transform は \(3L-3F\) である。

\(\delta=1\) の横断型は node 又は cusp のみである。
cusp の導手は \(2\Gamma\) なので、\(3L-3F\) の奇数係数に矛盾する。
従って全残余で generic node となり、導手は reduced。
generic node の二 branch は、大域的には一つの prime（二重 cover も含む）か、二つの prime である。
特殊 codimension 二の導手を消す議論ではなく、全 divisorial conductor の議論である。

## 2. 閉留数一形式：single prime の排除と fixed D 上の一意性

被約の実退化因子を持つ symplectic bivector の逆二形式は \(D\) に simple pole を持ち、
その Poincaré residue \(\alpha\) は \(D_{\rm reg}\) 上閉一形式である。
node の normalization branch へ Poisson bracket は regular に持ち上がり、導手に tangent である。
これは branch の minimal prime ideal を Hamiltonian derivation が保存することから従い、
係数体の元が Poisson central であることは仮定しない。
従って \(X\) の全 codimension 一で

\[
 \alpha=\iota_{\pi_X}\Omega_X
\]

は高々導手に simple pole を持つ。
導手の一般点で \(\Omega_X=(du/u)\wedge ds\wedge dt\) と書けば、tangency により
\(\alpha=c\,du/u+\text{regular}\) となる。
\(\pi|_S\) の一般 rank が二なので \(c\ne0\)。閉性から各 prime 上の residue は非零定数である。

導手が single prime \(\Gamma\) なら \(L^2\Gamma=12\)。
一般二 hyperplane の pullback curve を、全 bad codimension 二 locus と vertex を避けて選ぶ。
その smooth compact curve 上の全 residue 和は \(12\lambda\) となるので、
residue theorem により \(\lambda=0\)、上の非零性に矛盾する。
connected non-split double cover の場合も同じ一つの prime 上の定数を数えるので排除される。

split の場合は \(C_X=\Gamma_1+\Gamma_2\)、各 prime が \(S\) へ degree 一で写り、
\(L^2\Gamma_i=6\)。全 residue は
\(6\lambda_1+6\lambda_2=0\) を満たす。
ここから fixed \(D\) 上の \(\pi\) の projective choice は高々零次元と分かる。
実際、同じ pole set について residue が零の閉一形式は \(X_{\rm reg}\) で holomorphic。
small resolution の exceptional curve を除いた部分へ pullback し、
smooth \(\widetilde X\) 上の locally free \(\Omega^1\) の Hartogs 延長で全体へ延びる。
\(\widetilde X\) は \(\mathbb P^1\) 上の rational \(\mathbb P^2\)-bundle だから
\(H^0(\widetilde X,\Omega^1)=0\)。
従って residue pair は閉一形式を単射的に決定し、その線形空間は高々一次元である。

\(\omega_D\simeq O_D\) の generator を固定すると、\(\pi_D\mapsto\iota_{\pi_D}\Omega_D\) は
generic smooth locus で単射である。さらに ambient restriction も単射である。なぜなら

\[
 H^0(\mathbb P^4,\Lambda^2T(-5))=H^0(\mathbb P^4,\Omega^2)=0
\]

だからである。従って同じ \(D\) を持つ全 ambient Poisson tensors の projective fiber は高々一点。
local smoothness、scheme の reducedness、全 tensor の orbit 一意性は結論していない。

## 3. 型 014 と 005 の全補集合

014 で split なら \(\Gamma_i\sim a_iL+b_iF\) とする。
fiber \(\mathbb P^2\) への restriction より \(a_i\ge0\)。\(a_i=0\) の prime は fiber だが次数一なので不可能。
従って \(a_i\ge1\)、\(a_1+a_2=3\)、\(5a_i+b_i=6\) より、交換を除いて

\[
 \Gamma_1\sim L+F,\qquad\Gamma_2\sim2L-4F.
\]

\(\widetilde X=\mathbb P(O\oplus O(1)\oplus O(4))\) の fiber 座標を \(u,v,w\) とすると、
\(\operatorname{Sym}^2(O\oplus O(1)\oplus O(4))(-4)\) の全 section は

\[
 w\,(a u+\ell_1(s,t)v+q_4(s,t)w).
\]

従って \(|2L-4F|\) に prime divisor はない。
両因子の divisor は非零で、small resolution は divisor を縮めないので、画像上でも prime にはならない。
single prime は §2 で排除済みであり、型 014 は空である。

005 の vertex は line であり、finite 線形射影はその line を同次元の vertex に写す。
像 \(D\) は plane quintic curve とこの line の join である。
vertex を除いた局所環は base curve と affine space の積だから、二次元非正規支持の各成分は
base curve の singular point と vertex line の join、すなわち plane である。
vertex 自体は一次元にすぎない。従って integral degree 六 \(S\) は存在しない。

## 4. 型 023 の全射影空間と導手の閉条件

\(X=X_{023}\) の vertex 座標を \(z_0\)、base scroll を
\(Y=S(2,3)\subset\mathbb P^6\) とする。
\(H^0(X,L)=\mathbb Cz_0\oplus H^0(Y,H)\) は八次元。
basepoint-free な五次元 subspace 全体は
\(U\subset\operatorname{Gr}(5,8)\)、irreducible open、\(\dim U=15\) である。
全 \(V\in U\) は vertex で非零なので

\[
 W=V\cap H^0(Y,H),\qquad\dim W=4.
\]

\(W\) は \(Y\) 上 basepoint-free である。
もし \(y\) で全 \(W\) が零なら、その coneline 上で残る section \(z_0+a\) も零にでき、
\(V\) の basepoint が生じるからである。逆向きも成立する。
従って \(U\to U_4\subset\operatorname{Gr}(4,7)\) は
\(H^0(Y,H)/W\) を fiber とする三次元 affine bundle、\(\dim U_4=12\)。
vertex translation \(z_0\mapsto z_0+a\) を使えば、像 \(D\) は
\(\nu_W:Y\to T_W\subset\mathbb P^3\) の像 quintic を base とする cone である。
この全 \(\nu_W\) は、ample \(H\) と素数次数五および四 section の独立性により有限双有理である。
vertex 外の polynomial extension は normalization と conductor を保存するので、
\(X\) の divisorial conductor の splitting は、\(Y\) 上の conductor の splitting と同値。
vertex 点から追加の divisor は生じない。

ここで全 conductor の \(W\) に関する algebraic variation を説明する。
universal 四 section の五次乗の restriction map の kernel は各 fiber で一次元である。
各像が quintic hypersurface で、五次未満に関係式を持たないためである。
従って kernel は line bundle、universal image は相対 quintic hypersurface \(\mathcal T\)。
\(Y\times U_4\to\mathcal T\) は proper かつ全 fiber が有限なので finite である。
source と target はともに相対 Gorenstein で、有限双対性から conductor line は

\[
 \omega_{Y\times U_4/U_4}\otimes\nu^*\omega_{\mathcal T/U_4}^{-1}
 \simeq O_Y(-3H+3f)\otimes\text{(base line)}
\]

となる。この双対性は各 fiber の有限双対性と一致する。
conductor ideal の inclusion から、局所的に base line を trivialize すれば
全導手を表す非零 section の algebraic family、すなわち
\(U_4\to\mathbb P H^0(Y,O_Y(3H-3f))\) を得る。

\(Y=\mathbb F_1\) で \(3H-3f=3C_0+6f\)。
nonintegral divisor の集合はこの complete linear system で閉である。
具体的には有効錐が \(C_0,f\) によって生成されるため、二つの非零 effective class への分解は
\(0\le a\le3\)、\(0\le b\le6\) の有限個だけ。
各分解の section の multiplication map の projective image は proper image であり閉である。
その有限和が reducible 又は nonreduced divisor の全体を与える。

`SCROLL23_CONDUCTOR_WITNESS.md` の一つの \(W\) では、**全** conductor が C 上 integral と証明された。
従って nonintegral、特に split conductor locus は \(U_4\) の proper closed subset で、次元高々 11。
その \(U\) への pullback は次元高々 14。
この properness は Poisson に必要な conductor に限った標本数えではなく、全有限射影空間の閉条件である。

## 5. Aut の自由 fiber を引いた全 G 次元

embedded \(X_{023}\) を保存する群には、次元 14 の連結部分群がある。
\(E=O\oplus O(2)\oplus O(3)\) の bundle automorphism の次元は

\[
 \dim\operatorname{End}(E)=3+3+4+2=12.
\]

scalar を除き、base \(\operatorname{PGL}_2\) の三次元を加えて \(12-1+3=14\)。
作用は small resolution の dense open で faithful なので、\(X\) 上でも faithful。
これは Aut 全体の上界を仮定する数え方ではなく、実在する部分群の下界で十分である。

五 section の subspace に target projective frame を付けた projection-map 空間は
\(\operatorname{PGL}_5\)-torsor、frame の次元は 24。
split conductor の map 空間は次元高々 \(14+24=38\)。
上の 14 次元群による precomposition は image \(D\) を変えない。
しかもこの作用は自由である。もし \(\nu\circ g=\nu\) なら、有限双有理性から function field 上で \(g=1\)、
従って \(g=1\) だからである。
従って全 split image divisors の次元は高々

\[
 38-14=24.
\]

single prime は §2 で存在せず、split における fixed \(D\) の projective Poisson fiber は高々零次元。
従って全型 023 の当該 Poisson 像も次元高々 24 である。
型 014、005 は §3 で空なので、命題が従う。

## 6. 仮定・補集合・独立監査の状態

全 node/cusp、single prime（connected non-split を含む）/split、singular scroll 三型、
全 basepoint-free 五 section と target frames を保持した。
\(S\) の normality、smoothness、全域 ordinary-node、clean crossing、vertical divisor 不在は仮定しない。
closedness と一 witness が与えるのは G の全像次元上界であり、split 全点の具体的正常形ではない。
旧 rank 零の cone 排除を rank 二へ拡張せず、新 GB・射影点列・高次 jet 反復も使っていない。
smooth normalization 二型の別証明と本命題の独立監査を統合するのは integrator の仕事である。
