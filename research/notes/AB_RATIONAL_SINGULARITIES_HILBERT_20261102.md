# ABの有理特異点と全Hilbert級数

2026年11月2日。被約ABの全域正規性を入力とする。
以下は次数を増やした実験による予測ではなく、局所分解から得る有限次数上界を使う証明である。

## 1. incidenceの特異集合をblowupする

既存incidence f:X→ABはZ=LG(2,4)上局所的にY×A6である。
Y={N∈sl3:rank N≤1}は次元4の標準次数錐で、その射影化は

    F={([u],[v])∈P²×(P²)* : v(u)=0}.

Fは滑らかな三次元の(1,1)超曲面で、Yの標準埋込みのO_F(1)はO_F(1,1)である。
Xの被約特異集合Cは局所的に{0}×A6である。
Cでのblowup g:Xhat→Xは局所的にBl_0(Y)×A6であり、
Bl_0(Y)=Tot O_F(−1,−1)となる。従ってXhatは滑らかで、gはprojective birationalである。
中心はXの内在的な被約特異集合なので、chart間でこの構成は貼り合う。

F上の正の全次数について、P²×(P²)*の完全列

    0→O(n−1,n−1)→O(n,n)→O_F(n,n)→0

を使う。n≥1では左右のambient line bundlesの高次cohomologyは全て零。
n=0でもO(−1,−1)はacyclic、Oの正cohomologyは零なので、同じ結論となる。
従って全n≥0でH^i(F,O_F(n,n))=0 (i>0)である。
line bundleの全空間のaffine射影を使うとR^ig_*O_Xhat=0 (i>0)が局所的に従う。
g_*O_Xhat=O_XはXの正規性による。

前ノートでRf_*O_X=O_ABを証明済みなので、合成の導来順像から

    R(fg)_*O_Xhat=O_AB.

fgは滑らかな始域からのproper birational射である。
従ってABは原点を含め有理特異点を持ち、特にCohen–Macaulayである。
有理特異点からCMへの標準的帰結は、一次文献
Karen E. Smith, *F-rational rings have rational singularities*,
Amer. J. Math. 119 (1997), 159–180、§1 Definition 1.6直後で確認した。
[著者公開本文](https://sites.lsa.umich.edu/kesmith/wp-content/uploads/sites/1309/2024/07/finalration.pdf)
を2026年11月2日に確認。ここで同論文の正標数の主定理を使うわけではない。

## 2. 全次数の分子を有限個の値から決められる理由

B=C[AB]、S=Sym(W45*)とする。前ノートの最小自由複体Fは、
正cohomology消滅によりBの最小graded自由分解となっている。
局所Torの全shiftは、Yの

    (homological degree, internal degree, rank)
    (0,0,1),(1,2,9),(2,3,16),(3,4,9),(4,6,1)

と34個の一次Koszul因子の畳み込みである。
従って全Torの内部次数は40以下、また内部次数−homological degreeは2以下である。
同じhypercohomology spectral sequenceの各subquotientで内部次数は保たれるため、
Bの自由分解の全shiftも40以下となる。

従ってHilbert級数をH_B(z)=Q(z)/(1−z)^45と書くとdeg Q≤40である。
dim B=10なのでQは(1−z)^35で割り切れる。
よって

    H_B(z)=h(z)/(1−z)^10,    deg h≤5.

この分子次数上界はCM性を使わず、全域正規性と局所Torだけから従う。
従ってBの次数0から5の値を計算すれば全Hilbert級数が決まる。
前ノートにより全H^i(Z,E_n) (i>0)が零なので、各B_nはEuler特性で計算できる。
既存core/vertical Levi weight、Yの斉次ideal、およびBott操作を用いると

    dim B_n (n=0,...,5) = 1,45,761,6691,38570,166818.

これらと(1−z)^10を掛けて

    h(z)=1+35z+356z²+986z³+715z⁴+91z⁵

を得る。追加でn=6,7,8のEuler値とこの級数を照合するが、これらの追加値は
全次数を正当化する根拠ではない。

## 3. 帰結と有限証明書

ABの原点のHilbert–Samuel multiplicity、およびProj Bの次数は

    h(1)=2184

である。Bは標準次数CM環なのでreg B=deg h=5、a(B)=5−10=−5、
S上のprojective dimensionは45−10=35となる。
またF^(−p)に寄与するTor_iはi=p+j、j≤3であり、
shift≤i+2≤p+5である。特にABの斉次定義idealは次数6以下で生成される。
次数二の関係の次元はbinom(46,2)−761=274である。
次数三から六の最小生成元数や、明示的な全定義idealを求めたという主張はしない。

有限証明書は全Koszul畳み込み、全shift上界、低次Euler値、分子と余分な三次数の照合、
Fのsection空間の寸法恒等式

    binom(n+2,2)²−binom(n+1,2)²=(n+1)³

を確認する。blowupの幾何とcohomology消滅の論証は本文にあり、プログラムだけの証明ではない。

この結果は被約ABの特異点を解決する。T01の実際の零fiberの厚み、
全S61の成分尽くし、一般アフィン五点以外の不存在には別の議論が必要である。

    .venv/bin/python code/ab_rational_singularities_hilbert_20261102.py --output reports/AB_RATIONAL_SINGULARITIES_HILBERT_20261102.json
