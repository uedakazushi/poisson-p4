# 底を優越する未知成分の次元と正規退化因子

2026-11-07。Gの次の幾何的入口。ここでは必要条件を証明し、次に必要な
非優越命題を一つ特定する。その非優越命題は未証明であり、Gの完了を主張しない。
交換子・多項式左逆・有限射影だけでは孤立残余を排除できないことは
[独立reviewの反例](NORMAL_DIRECTION_INCIDENCE_REVIEW.md)による。

## 1. 未知成分の射影次元は25、26、27のいずれか

\(W\simeq\mathbb A^{16}\)をprimitive一次部分の空間、
\(S\to W\)を元の全134式によるS61とする。以下の次元はPGL商の前である。
既知五点を分離した一般ファイバーの残余代数を\(B\)と書く。

**命題1。** \(B_{\mathrm{red}}\)が空でなければ、それに対応する全Poisson空間の
未知既約成分\(C\)が存在し、ある\(d\in\{0,1,2\}\)について

\[
 \dim C=25+d
\]

となる。ここで\(d\)は対応する一般\(W\)-ファイバーの次元である。

**証明。** 一般ファイバーの既約成分の閉包は\(W\)を優越するS61の成分を与える。
既存の[アフィン次元上界](../notes/S61_AFFINE_SURFACE_BOUND_20261102.md)により
\(d\leq2\)であり、この成分の次元は\(16+d\)である。
[C0の大域gauge監査](C0_S61_SLICE_SCALAR_AUDIT_20261106.md)の第1節を用いる。
固定点における非退化な射影評価値は5次元、primitive gaugeで除いたperspectiveは
4次元なので、全射影Poisson空間の対応成分には合計9次元が加わる。
評価値の安定化群CSp4とperspective群は連結なので、有限個の成分を置換しない。
従って\(\dim C=(16+d)+9=25+d\)である。
これはS61の55消去とBianchiの大域恒等式を入力とするC0と同じ依存を持つ。□

[既知一次像の監査](EXCEPTIONAL_BASE.md)によれば、既知の成分で\(W\)を優越するのは
T01、T39、T40のみで、その一般ファイバーが既知五点である。
従って未知\(C\)はこれらの再表示ではない。

MPSのTheorem 1.4はSNC log-symplectic locusと交わる成分を40個に尽くす。
したがって未知\(C\)はSNC log-symplectic点を一つも含まない。
一般元だけが非SNC、という弱い主張ではない。
一次資料は[MPS, arXiv:2010.08692v2](https://arxiv.org/html/2010.08692v2)、
Theorem 1.4、Example 6.24。既知40との対応は
[成分監査](KNOWN_COMPONENT_AUDIT.md)第4節を用いる。

## 2. 正規被約退化因子ならunimodular零点が必要

\(\mathcal P\subset\mathbb P H^0(\mathbb P^4,\Lambda^2T)\)をPoisson空間とし、
\(D_\pi=V(\pi\wedge\pi)\)を退化五次因子とする。
\(\mathcal N\)を一般階数4かつ\(D_\pi\)が被約正規である部分とする。
被約Pfaffianを持つPoisson構造はlog-symplecticである。

局所体積形式\(\mu\)に対するmodular vector fieldを\(Z_\mu\)と書く。
\(\pi(p)=0\)なら\(Z_\mu(p)\)は体積形式に依らない。
実際、体積形式の変更で加わるHamiltonian vector fieldは\(p\)で零になる。
この値を\(Z_\pi(p)\)と書く。

**命題2。** 未知成分\(C\)の一般退化因子が被約正規なら、その一般元\(\pi\)は

\[
 \pi(p)=0,\qquad Z_\pi(p)=0
 \tag{1}
\]

を満たす点\(p\in\mathbb P^4\)を持つ。

**証明。** 正規log-symplectic因子の特異 locusは純余次元3のPoisson部分空間である。
4次元ではその各点を通るsymplectic leafは零次元なので、そこでは\(\pi=0\)である。
全ての特異点で\(Z_\pi\ne0\)なら、PymのDefinition 4.2とRemark 4.4により
\(\pi\)はpurely ellipticである。PymのTheorem 1.3により、この構造は
\(q_{5,1}\)、すなわちT39に属する。未知成分の一般元にはこれは起こらない。
従ってある特異点で(1)が成り立つ。□

使用する一次資料は[Pym, *Elliptic singularities on log symplectic manifolds
and Feigin–Odesskii Poisson brackets*](https://arxiv.org/pdf/1507.05668)、
Theorem 3.1、§4.1、Definition 4.2、Remark 4.4、Theorem 1.3である。
今回これらの適用仮定を本文で照合した。正規性だけでpurely ellipticとはしない。
同論文が明示する通り、\(q_{5,2}\)の退化因子は正規でないため、この分岐の対象外である。

## 3. 必要条件をproperな閉incidenceとして定式化する

\[
 \mathcal I_{\rm um}
 =\{([\pi],p)\in\mathcal P\times\mathbb P^4:
        \pi(p)=0,\ Z_\pi(p)=0\},
 \qquad
 \mathcal U=\operatorname{pr}_{\mathcal P}(\mathcal I_{\rm um})
 \tag{2}
\]

と定義する。\(\mathcal I_{\rm um}\)は閉であり、射影がproperなので
\(\mathcal U\subset\mathcal P\)も閉である。
閉性は任意のaffine座標で\(\pi(p)=0\)を課した後、一次jetの自然な縮約

\[
 T_p^*\otimes\Lambda^2T_p\longrightarrow T_p
\]

が零であるという線形条件から従う。これは局所体積形式や\(\pi\)の非零スカラー倍に
依らず、全射影空間で貼り合う。座標表示では、符号の規約を除いて
\(Z^j(p)=\sum_i\partial_i\pi^{ij}(p)\)である。

命題2と閉性から、一般退化因子が正規被約である未知成分は全体として
\(C\subset\mathcal U\)を満たす。
ただし\(\mathcal U\)は未知成分専用の集合ではなく、既知成分の境界も含み得る。
同じ証明は点ごとの包含
\(\mathcal N\setminus C_{T39}\subset\mathcal U\)も与える。

さらに(1)には以下の有限jetの必要条件がある。

**命題3。** (1)の点を原点とし、\(\pi=\pi_1+\pi_2+\cdots\)をTaylor展開する。
\(\pi_1\)は4次元unimodular Lie algebraの線形Poisson構造であり、
\(\pi_1\wedge\pi_1=0\)である。従って\(\operatorname{mult}_p D_\pi\geq3\)。
局所Pfaffianの三次斉次項を\(f_3\)とすれば
\(\pi_1^\#(df_3)=0\)である。

**証明。** Jacobiの一次斉次部分から\(\pi_1\)は線形Poissonであり、
その構造定数が\(T_p^*\)上にLie bracketを定める。
\(Z_\pi(p)=0\)は\(\operatorname{tr}(\operatorname{ad}a)=0\)を全\(a\)について
課す条件、すなわちunimodularityである。
\(f_2=\operatorname{Pf}(\pi_1)\)とする。
Poisson構造のmodular恒等式\(\pi_1^\#df_2=f_2Z_{\pi_1}\)で右辺は零である。
\(f_2\ne0\)なら\(D(f_2)\)上で\(\pi_1\)が可逆なので\(df_2=0\)。
標数零で\(f_2\)は正次数斉次式だから矛盾する。
従って\(f_2=0\)で、退化因子の次数2以下の項は消える。
全\(\pi\)の同じmodular恒等式を次数3で比較すると、
\(Z_\pi\)の定数項が零であるため\(\pi_1^\#df_3=0\)を得る。□

このLie algebra条件は必要条件である。
Lie bracketの全型を列挙しただけでは、その全高次延長や全Poisson構造を被覆したことにはならない。

## 4. 非可解一次Lie代数の分岐は全高次延長ごと排除できる

**命題4。** 滑らかな複素4次元多様体上のPoisson構造が\(p\)で零となり、
その一次Lie代数が\(\mathfrak{sl}_2\oplus\mathbb C\)なら、
\(p\)の近傍で\(\pi\wedge\pi=0\)である。
従って一般階数4のPoisson構造の零点では、一次Lie代数は必ず可解である。
この主張には退化因子の正規性やmodular vectorの消失を仮定しない。

**証明。** \(\mathfrak g=\mathfrak{sl}_2\oplus\mathfrak z\)、
\(\dim\mathfrak z=1\)とする。一次Poisson構造\(\pi_1\)に対する、係数が
\(k\)次斉次のpolyvectorの複体は

\[
 C^\bullet(\mathfrak g,\operatorname{Sym}^k\mathfrak g)
\]

というChevalley–Eilenberg複体である。
中心\(\mathfrak z\)は\(M_k=\operatorname{Sym}^k\mathfrak g\)に自明に作用するので、
複体自体が\(C^\bullet(\mathfrak{sl}_2,M_k)\otimes\Lambda^\bullet\mathfrak z^*\)
と分解する。有限次元加群\(M_k\)にWhiteheadの二補題を用いると

\[
 H^2(\mathfrak g,M_k)
 =H^2(\mathfrak{sl}_2,M_k)
  \oplus H^1(\mathfrak{sl}_2,M_k)\otimes\mathfrak z^*=0.
 \tag{3}
\]

中心が1次元で\(\Lambda^2\mathfrak z^*=0\)であることがここで必要である。

形式座標で\(\pi=\pi_1+\pi_k+\text{次数}\geq k+1\)、\(k\geq2\)まで
低次項を消したとする。Jacobiの次数\(k\)部分は
\([\pi_1,\pi_k]=0\)である。(3)より係数次数\(k\)のvector field \(X_k\)があり、
\(\pi_k=[\pi_1,X_k]\)とできる。\(X_k\)の形式flowを適切な符号で施せば
\(\pi_k\)は消え、それより低い次数を変えない。
これを全\(k\)について行う合成は、接写像が恒等である\(\mathfrak m\)-進自己同型を定める。
従って任意の全高次Poisson延長は形式的に\(\pi_1\)へ線形化される。

\(\pi_1\)は\(\mathfrak{sl}_2^*\)の3方向だけに支えられるので
\(\pi_1\wedge\pi_1=0\)。形式座標変換を戻すと、元の\(\pi\wedge\pi\)も完成局所環で零である。
滑らかな代数的または解析的局所環からその完成への写像は単射なので、実際のgermでも零となる。
この段階で形式座標変換の解析的収束は不要である。
特に既約な\(\mathbb P^4\)上では\(\pi\wedge\pi\)は全域で零となり、一般階数4に反する。

最後に、4次元の非可解複素Lie代数のLevi因子は\(\mathfrak{sl}_2\)で、radicalは1次元である。
\(\mathfrak{sl}_2\)からその1次元radicalの導分代数への作用は零なので、
Lie代数は\(\mathfrak{sl}_2\oplus\mathbb C\)である。これで非可解の場合を全て排除した。□

形式線形化は既知のLevi理論の帰結であり、新規性を主張しない。
一次資料の[Zung, *Levi decomposition of analytic Poisson structures and Lie algebroids*,
arXiv:math/0203023v2](https://arxiv.org/pdf/math/0203023)、Theorem 1.1直後のRemark 3は、
半単純代数に1次元中心を加えた場合の解析的線形化も明記している。
今回その本文を確認した。上の証明では分類に足りる形式的結論だけをWhiteheadから導いた。

この排除は任意の零点での内在的なLie代数を扱うため、S61のrank4基点を変えても成立する。
Poisson高次項を切り捨てず全次数を扱っており、一つの係数切片だけの排除ではない。
一方、可解unimodular一次部の全Poisson延長は未処理である。

## 5. 次に必要な一つの命題と残る補集合

S61の標準化された構造を\(\pi_s\)とし、

\[
 J_{\rm um,nor}=
 \{s\in S: [\pi_s]\in(\mathcal N\cap\mathcal U)\setminus C_{T39}\}
\]

と置く。既知T39閉包を除くことで、その境界の解析を余計な義務にしない。
次の**未証明命題**がこの正規被約分岐を閉じるための具体的な標的である。

\[
 \dim\overline{\operatorname{im}(J_{\rm um,nor}\longrightarrow W)}\leq15.
 \tag{G-normal}
\]

これが証明できれば、十分一般の\(\beta\)における正規被約退化因子を持つ解は
全てT39に属する。証明は上の点ごとの包含であり、未知fiberの次元や非被約構造を
さらに決める必要はない。
より強い全空間上界\(\dim(\mathcal U\cap\mathcal N)\leq24\)も命題1により十分だが、
底像の上界で足りるので、強い方を前提に課さない。

(G-normal)の証明候補は(2)を零点\(p\)における**可解unimodular**一次Lie bracketで層化し、
全Poisson延長を保持したまま各層のprimitive一次像を評価することである。
rank4の基点で作る\(\beta\)と、特異点\(p\)で作るLie bracketは異なるjetである。
基点を動かす操作と全高次条件を省略して、後者の型の個数から前者の次元を推定してはならない。
このnoteにはその像上界の証明や計算結果はまだない。

また(G-normal)が閉じても、一般退化因子が**非被約**の未知成分と、
**被約だが非正規**の未知成分は残る。T40自体が後者に入る。
この二つの補集合や特殊\(\beta\)の問題を、正規分岐の議論で排除したとはしない。

## 独立確認の範囲

2026-11-07、別担当が命題1のC0次元橋、命題3のmodular恒等式と次数比較を読解確認した。
Pym一次PDFのDefinition 4.2、Remark 4.4、Theorem 1.3を別途読み、命題2の適用を確認した。
命題4についても別担当が複体の直和分解、中心1次元の必要性、全次数の形式帰納、
完成からの零性の帰還、および4次元の非可解Lie代数の帰着を読解確認した。
全134式や旧境界証明書の再計算は行っていない。(G-normal)はこの確認後も未証明である。
