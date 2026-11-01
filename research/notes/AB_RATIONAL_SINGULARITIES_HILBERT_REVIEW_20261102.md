# AB の有理特異点と Hilbert 級数：独立レビュー

2026年11月2日。[主ノート](AB_RATIONAL_SINGULARITIES_HILBERT_20261102.md)を受理する。
被約 AB は有理特異点を持ち、Cohen–Macaulay である。その Hilbert 級数は

    (1+35z+356z²+986z³+715z⁴+91z⁵)/(1−z)^10

となる。以下は既存の全域正規性を依存とする証明監査と、別構成の有限係数検算である。

## 1. 解消と消滅

incidence X は LG の chart 上で Y×A³×A³ であり、
Y は trace 零 rank≤1 の 3×3 行列錐である。特異点集合は Y の原点に対応する
滑らかな閉集合 C。C の blowup は、各 chart で Bl_0(Y)×A⁶ となる。
Y は標準次数環の錐なので、その原点 blowup は、滑らかな incidence
F⊂P²×(P²)* 上の O_F(−1,−1) の全空間である。
この局所記述は内在的な中心 C の blowup の chart であり、独立に選んだ解消を
貼り合わせるという仮定を要しない。従って Xhat は滑らかである。

F の (1,1) 超曲面完全列を n≥0 で使う。n=0 の左端 O(−1,−1) も acyclic なので、
H^i(F,O_F(n,n))=0 (i>0) は例外なく成立する。
line bundle 全空間の affine 射影で O の順像は ⊕_{n≥0}O_F(n,n)。
F は射影的なので有限 affine cover による Čech 計算で cohomology とこの直和が可換であり、
blowup の高次順像は零になる。零次順像は X の正規性から O_X である。
これと既証明 Rf_*O_X=O_AB を合成して、滑らかな proper birational 解消の
正の高次順像が全て零と分かる。

有理特異点から Cohen–Macaulay 性への帰結は主ノートの引用先を実際に確認した。
Karen E. Smith, *F-rational rings have rational singularities*, §1 の有理特異点の
定義直後にこの性質が明記されている。
[著者公開本文](https://sites.lsa.umich.edu/kesmith/wp-content/uploads/sites/1309/2024/07/finalration.pdf)。
正標数の主定理を用いる議論ではない。

## 2. 全次数の上界

既存の全域正規性レビューで用いた最小 graded perfect complex は、いま B の最小自由分解である。
その導来零 fiber は H^j(LG,Tor_i) の spectral sequence の subquotient で計算される。
局所 Tor の内部次数は、Y の shifts (0,2,3,4,6) と34個の一次 Koszul 因子を畳み込んで
高々40である。内部次数は全ての differential で保存されるため、global な分解の shift も高々40。
Tor 束の大域分裂や spectral sequence の退化は不要である。

従って H_B=Q/(1−z)^45 において deg Q≤40。
Krull 次元10から z=1 における零の位数が35なので、Q は (1−z)^35 で割り切れ、
既約分母 (1−z)^10 に対する分子の次数は高々5となる。
この段階では CM 性を用いていない。従って有限個の Hilbert 値から全級数を決定する論法は妥当。

さらに total degree −p に寄与する局所 Tor は i=p+j、0≤j≤3 であり、
d−i≤2 から d≤p+5。p=1 の関係式は次数6以下で生成される。
CM 性と dim B=10 を合わせると pd_S B=35、reg B=5、a(B)=−5 を得る。
次数2184は分子の値 h(1) である。

## 3. 独立 Cauchy 計算

主 producer の Gröbner 単項式 character や Bott 関数を import していない。
rank 一行列の Segre 表示から、E=Sym²(C²) として core の次数 n character を

    det^n [Sym^n(E)⊗Sym^n(E*) − Sym^(n−1)(E)⊗Sym^(n−1)(E*)]

で直接組み立てる。第二項は trace の一式を除いたもの。
これを vertical weights (3,1),(2,2),(1,3) の対称冪と畳み込む。
GL2 最高重みの重複度は隣接する weight multiplicity の差から求め、
Euler 特性は符号を含む Weyl 多項式

    χ(a,b)=(a−b+1)(b+1)(a+2)(a+b+3)/6

を使う。この方法は主の dominant chamber への移動実装と異なる。
再構成した全 weight character が元の character に戻ることも確認する。

独立結果は n=0,…,8 について

    1,45,761,6691,38570,166818,586770,1766622,4711707

であり、とくに未使用だった次数5の値166818を別計算で得た。
これらは associated graded bundle の Euler 特性である。
実際の B_n の寸法への同定には、既証明の正 cohomology 消滅を用いる。
次数0,…,5と分子次数上界だけで結論は決まり、6,…,8は追加照合にすぎない。

    .venv/bin/python code/ab_rational_singularities_hilbert_review_20261102.py --output reports/AB_RATIONAL_SINGULARITIES_HILBERT_REVIEW_20261102.json

Y の自由分解 exactness、既存 incidence の同型、全域正規性と高次順像消滅は凍結済みの依存である。
本レビューはそれら全てを別実装したという意味ではない。
T01 のスキーム論的零 fiber の厚み、S61 の尽くし、成分数の上界は未解決のままである。
