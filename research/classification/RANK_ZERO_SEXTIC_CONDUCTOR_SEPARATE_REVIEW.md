# 次数六 rank 零曲面の全排除：別作者による導手・Poisson 持上げの監査

2026-11-07。producer は `generic_exhaustion`、reviewer は `log_manuscript`。
主義務 G の有限な幾何証明監査であり、新しい CAS 計算ではない。
`RANK_ZERO_SEXTIC_CONDUCTOR_EXCLUSION.md` の結論を、以下の全仮定・全補集合の下で受理する。
対応する構造化記録は
`reports/CLASSIFICATION_RANK_ZERO_SEXTIC_CONDUCTOR_REVIEW_20261107.json`、
監査契約は `jobs/rank_zero_sextic_conductor_review.json` である。

**受理する命題。** 複素射影四次元空間上の rank 四 Poisson tensor の実退化因子
\(D\) が既約被約 quintic であるとする。その非正規支持の integral surface 成分
\(S\) が次数六で、\(\pi|_S=0\) となることはない。
\(S\) の normality、smoothness、cone 性、一般 ordinary-double 型、
あるいは vertical divisor の不存在は仮定しない。

## 1. 一次資料の範囲と既受理入力

Andreas Höring, *The sectional genus of quasi-polarised varieties*,
[author PDF](https://math.univ-cotedazur.fr/~hoering/articles/a12-sect-genus.pdf)
の Theorem 1.2 とその直後、印刷二頁を reviewer が直接読解した。
同定理の仮定は normal complex projective variety、nef and big Cartier divisor、
sectional genus 零である。smooth、terminal、klt 又は Cohen--Macaulay は追加されていない。
従うのは \(\Delta=0\)、すなわち \(h^0(L)=\dim X+L^{\dim X}\) である。
直後の段落は complete linear system による normal な very ample model も記す。
本件では \(L\) が ample なので、この model への射は finite でもあり同型となる。

minimal-degree variety の全分類は Brodmann--Park,
[arXiv:0808.0090v3](https://arxiv.org/pdf/0808.0090v3),
Notation and Remark 2.1、印刷四頁で別に直接照合した。
ここでは同論文の outer *point* projection のみの定理を、
本件の任意の三次元中心の射影へ拡大していない。
使用するのは minimal-degree 分類そのものである。

既受理の数値入口は `IRREDUCIBLE_NONNORMAL_QUINTIC_ENTRANCE.md` §1、
cone の全 Euler 源は `RANK_ZERO_HIGH_DEGREE_CONE_EXCLUSION.md`、
退化曲面の prime-CI-square 排除は
`IRREDUCIBLE_NONNORMAL_DOUBLE_SURFACE_DEGREE_RESTRICTIONS.md` §1 である。
これらの適用仮定も読み直した。旧 normal-sextic の ruled 数値式や未読分類は使わない。

## 2. 二つの Bertini 切断と sectional genus の一致

正規化 \(\nu:X\to D\)、\(L=\nu^*O_D(1)\) とする。
\(X\) は normal threefold で、\(L\) は finite pullback だから ample、かつ元の五 section で生成される。
normality により \(\dim\operatorname{Sing}X\le1\) であり、\(\nu\) は finite なので
一般の平面 \(P\subset\mathbb P^4\) は \(\nu(\operatorname{Sing}X)\) を避ける。
二つの一般 hyperplane を引き戻す basepoint-free 系へ Bertini を順次適用すると、
\(\Gamma=\nu^{-1}(P)\) は全体が smooth locus 内の滑らかな曲線である。
二切断は regular sequence でもある。具体的には normal の \(S_2\) 性により
第一 Cartier 切断は \(S_1\)、第二一般切断はその有限個の associated generic point を避ける。

一般の平面 quintic \(C=D\cap P\) は integral である。
\(\Gamma\to C\) は finite で、\(C\) の一般点では \(\nu\) が同型だから birational である。
追加の曲線成分があれば finite 性によりそれも \(C\) 全体へ写り、一般点の一対一性に反する。
したがって \(\Gamma\) は integral な \(C\) の正規化である。
\(X\) 全体の CM 性も smoothness も、この議論には必要でない。

数値入口の
\[
 \sum_i\delta_i\deg S_i=6-g,\qquad L^3=5,\qquad K_X+C_X\sim0
\]
に \(\deg S=6\) を入れると、\(g=0\)、\(\delta_S=1\)、
二次元の非正規支持は \(S\) 一つだけとなる。
ここで \(g=p_a(\Gamma)\) は Höring の sectional genus と同じである。
二つの Cartier 完全列を Euler characteristic に適用すると
\(\Delta_t^2\chi(X,tL)=\chi(\Gamma,tL)=5t+1-g\) となるからである。
Höring の rising-binomial 表示では右辺の定数項が \(\chi_2(X,L)\) であり、
\(g(X,L)=1-\chi_2(X,L)=g\) である。
一般 plane curve の算術種数六と、この正規化種数零を混同していない。

## 3. 全 finite projection と五つの scroll

Höring の定理から \(h^0(X,L)=8\)。complete linear system は元の五 section を含み、
射 \(\varphi:X\to\mathbb P^7\) は \(\nu\) を linear projection として因数分解する。
\(L\) の ample 性から \(\varphi\) は finite、\(\nu\) の birational 性から \(\varphi\) も birational。
像は非退化で次数五、次元三、codimension 四の minimal-degree variety である。
上記 normal-model の結論、又は minimal-degree scroll の normality により、
\(\varphi\) は同型であり complete embedding となる。
よって元の \(\nu\) は、この \(X\subset\mathbb P^7\) の外に中心を持つ五 section の射影である。
中心を一つの点と仮定してはいない。

次数・次元から全型は
\[
 S(1,1,3),\ S(1,2,2),\ S(0,1,4),\ S(0,2,3),\ S(0,0,5)
\]
に尽くされる。quadric、projective space、Veronese とその cone は次数五でない。
正の三 summand の最初の二型は smooth、残る三型は vertex を持つ cone である。

cone 型の vertex \(V\) は \(X\) に含まれるため射影中心を避ける。
よって射影は \(V\) 上で単射で、その像 \(V'\) は同次元の linear vertex となり、
join の線形射影により \(D\) 自身も cone である。
cone の vertex を除いた局所 chart は base の affine chart と affine space の積である。
normality と normalization defect の支持はこの積で保たれるので、
唯一の二次元非正規支持 \(S\) は \(V'\) との join、すなわち cone である。
非退化なら高次数 cone の全源補題、hyperplane 内なら次数三以上の prime-CI-square 補題が使え、
いずれも実既約 quintic に反する。vertex が line となる最後の型も除外されず、
その場合の二次元 join は plane となることからも次数六に反する。

## 4. 導手公式と全 \(\delta=1\) 横断型

有限双対性を codimension 一で適用すると、Gorenstein quintic の
\(\omega_D\simeq O_D\) に対し
\(\omega_X\simeq\operatorname{Hom}_D(\nu_*O_X,O_D)\) を rank-one reflexive 層として得る。
これは conductor ideal であるため \(K_X+C_X\sim0\) となる。
singular \(X\) に Cartier canonical class を勝手に課してはいない。
cone 三型を除いた後は \(X\) が smooth scroll なので、ここから Cartier/Pic 計算が可能である。

smooth 二型で \(L\) と ruling fiber divisor \(F\) は Pic の自由基底で、
\(\det(O(a)\oplus O(b)\oplus O(c))=O(5)\) だから
\[
 K_X=-3L+(5-2)F=-3L+3F,\qquad C_X\sim3L-3F.
\]
finite map の divisor は曲面へ写るので、全 conductor prime はこの唯一の \(S\) 上にある。

\(S\) の一般点で局所化・完成し residue field を分離拡大すると、
\(D\) は reduced plane-curve singularity になる。
\(\delta=1\) なら multiplicity は二で、char 零 Weierstrass 正規形は
\(v^2-u^r\)、\(r=2\) 又は三である。
\(r\ge4\) なら \(\delta=\lfloor r/2\rfloor\ge2\) となる。
したがって全型は node 又は cusp であり、ordinary 性は仮定ではなく数値からの結論である。

cusp では \(O=K[[t^2,t^3]]\subset K[[t]]\) の conductor は \((t^2)\)。
よって \(C_X=2\Gamma\) であり、\(S\) 上の唯一の divisor の係数は二である。
\(X\) は smooth だから \(\Gamma\) は Cartier であるが、
Pic class \(3L-3F\) は二で割れない。codimension 二の特殊導手や vertical 部分を捨てなくても、
この divisorial parity の矛盾は変わらない。

node では split なら二 branch、non-split なら分離拡大後の二 branchであり、
全ての conductor valuation は一である。従って大域的に \(C_X\) は reduced divisor となる。

## 5. node での持上げ：係数体を中心と仮定しない直接証明

実退化因子の ideal は Hamiltonian derivation に保存される。
したがって \(D\) 上に bracket が降り、smooth dense open の bivector を
\(X\) の rational bivector \(\pi_X\) として考えられる。
正規化全体への一般延長定理を入力にはしない。

conductor の一般点で、finite separable residue extension と completion 後の環を
\[
 A=K[[u,v]]/(uv),\qquad \mathfrak m=(u,v)
\]
とする。\(K\) は \(S\) の function field の有限分離拡大であり、
\(K\) が Poisson 中心とは仮定しない。
\(\pi|_S=0\) は全ての bracket について \(\{A,A\}\subset\mathfrak m\) を意味する。
特に \(\mathfrak m\) は Poisson ideal である。
Leibniz 則から \(\{\mathfrak m^r,A\}\subset\mathfrak m^r\) となるので、
bracket は completion へ連続に延長する。
有限分離拡大の際も、primitive element の separable equation を微分し
unit であるその導関数で解くことで derivation と bracket は一意に延長し、
residue bracket 零性も保たれる。

任意の Hamiltonian derivation \(\partial\) について
\(v\partial u+u\partial v=0\)。これを \((u)\) で割った
\(K[[v]]\) では \(v\) が非零因子だから \(\partial u\in(u)\)。
同様に \(\partial v\in(v)\) であり、**各 minimal prime は個別に不変**である。
よって二つの全 branch 商
\(A/(v)=K[[u]]\)、\(A/(u)=K[[v]]\) に bracket が降りる。
どちらも環の全ての元を含む surjective 商であるため、
誘導 bracket の全値はそれぞれ \((u)\)、\((v)\) に入る。
これは tangential bracket だけ、又は選ばれた二座標だけの零性ではない。

より具体的には、\(K/\mathbb C\) の separating transcendence basis を二つ選び、
branch parameter と合わせると smooth normalization の generic divisor 付近の
微分基底となる。その三対の bracket は全て uniformizer で割れる。
従って \(\pi_X\) はそこで regular で、全 conductor branch 上で一次以上消える。
faithfully flat な分離拡大・completion でのこの valuation 条件は元の branch へ降下する。
連結二重被覆や branch exchange を仮定で排除していない。

\(C_X\) 以外の全 codimension 一では normalization は同型で、
normal な \(D\) はその codimension 一で regular だから \(\pi_X\) は regular である。
したがって smooth \(X\) 上で rational section を全 prime divisor で検査することにより
\[
 \pi_X\in H^0(X,\Lambda^2T_X(-C_X)).
\]
これは locally free 層の Hartogs 延長であり、singular surface \(S\) 上の
reflexive 形式を resolution へ延長する主張ではない。

三次元の恒等式と \(K_X+C_X\sim0\) から
\(\Lambda^2T_X(-C_X)\simeq\Omega_X^1\)。
smooth scroll \(p:X\to\mathbb P^1\) の相対 cotangent 完全列で、
\(p_*\Omega_{X/\mathbb P^1}^1=0\) は \(\mathbb P^2\) fiber 上の零性から従い、
\(H^0(\Omega_{\mathbb P^1}^1)=0\) と合わせて \(H^0(X,\Omega_X^1)=0\) となる。
一方 \(\pi_X\ne0\) である。実際、\(D\) の smooth 一般点で \(\pi\) は接する rank 二であり、
そこで全係数が零なら actual Pfaffian は \(D\) を二回含んで reduced 性に反する。
この非零 section と上の零空間が矛盾し、node も尽くされる。

## 6. 判定と残余

以上により、producer の命題を **全 degree 六 rank 零非正規支持の排除として受理**する。
旧 normal-sextic 排除の仮定を非normal 曲面へ広げたのではなく、
quintic 自身の normalization の sectional genus 零から別に証明した。
計算による rank、有限標本、formal jet、未証明の flatness 又は正規化上の一般 Poisson 延長を使わない。

この結論だけから degree 五、plane、quadric cone、一般 rank 二 restriction、
全 G 非優越、全 E 帰属又は分類定理は従わない。
既知成分数は変更しない。保存済み102成分主TeX/PDFは今回の監査では変更しない。
