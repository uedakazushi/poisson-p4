# degree四rank零曲面：全integral非退化支持の有限還元

2026-11-07、主義務G。本文はprojection packet間の網羅性を明示する。
`PROJECTED_QUARTIC_CONE_RANK_ZERO_EXCLUSION.md` の全source及びQ因子証明を使う。
そのproducer証明と本有限還元の独立review receiptは別noteで管理する。

**結論。** S⊂P4を任意のintegral・非退化degree四曲面とする。
rank四Poisson tensor πがS上で零となり、実anticanonical Dが既約被約quinticであると仮定する。
全支持は下の有限還元で覆われる。
SがCI22である場合の全sourceのPGL5-saturationは既証明のdimension≤24であり、
残る支持ではそのようなπは存在しない。
従ってこの分岐から未知のdimension≥25 componentは生じない。
normality、smoothness、ordinary-double型、vertical divisorの無いことを仮定しない。

## 1. 使用する外部分類と正確なscope

M. Brodmann–E. Park, *On varieties of almost minimal degree I: Secant loci of
rational normal scrolls*, arXiv:0808.0090v3,
<https://arxiv.org/pdf/0808.0090v3> を直接読解した。
Introductionの印刷1–2頁は、integral・非退化なalmost-minimal-degree variety全体を
normalかつlinearly normalなdel Pezzoの場合と、一つ上の射影空間にあるminimal-degree
varietyの外点からの単純射影の場合に分ける。
後者はnon-normal又はnon-linearly-normalの全場合を含む。
Remark6.1(A)、印刷16頁は、linearly normalなdel Pezzo varietyがACMであることを明記する。
minimal-degree分類はNotation and Remark2.1、印刷4頁で確認した。

これはgeneric centreだけの分類ではない。
surface degree四、ambient P4では、normal・linearly normalの場合はACMである。
それ以外の全Sは、degree四minimal-degree surface Y⊂P5を外点pから射影したimageである。
minimal-degree分類でYはVeronese surface、S(1,3)、S(2,2)、S(0,4)の四種類だけとなる。
dimension二なので、これより高次元のVeronese coneは入らない。

## 2. 全ACM支持はCI22へ入る

ACM degree四codimension二integral surfaceのcoordinate ringを三つのlinear parametersで割る。
Artin quotientはembedding dimension二、length四で、Hilbert functionは(1,2,1)である。
従ってI_S(2)は二次元。
独立な二quadricにlinear common factorがあれば、integralityによりSはhyperplaneに含まれ、
非退化性に反する。よって二quadricはregular sequenceでdegree四CIを定める。
このCIにdegree四Sが含まれるので、CM性・unmixednessとdegree比較によりschemeとして一致する。

全CI22は、pencilがnondegenerate memberを持つ場合と全memberがdegenerateな場合に分かれる。
前者の全Segre型は `REGULAR_PENCIL_CI22_ABELIAN_CONORMAL_NONDOMINANCE.md`、
後者の全Kronecker型は `ALL_DEGENERATE_CI22_ABELIAN_CONORMAL_EXCLUSION.md` で覆われる。
どちらも別作者review済みで、全rank零G像のdimension≤24を与える。
従ってnormal・linearly normal S及び後述のACM projection centresは全て既済枝へ入る。

## 3. Veroneseからの全centre

centreがVeronese surface自身にある場合は射影定義の外である。
closed secant variety内の外点ではimageはCI22である
（同一次資料Remark6.3(A)、印刷19頁）。
closed secant外では射影は同型埋込みで、smooth projected Veronese surfaceとなる。
その全rank零sourceからrank四Poisson tensorが存在しないことは既済のprojected-Veronese-zero
証明書・独立reviewが与える。

## 4. smooth quartic scrollからの全centre

`PROJECTED_QUARTIC_SCROLL_ISOLATED_MODELS.md` は次をwhole-centre orbitとして分類する。

| minimal-degree Y | 既済CI22へ入るcentre | 残る全centre orbit |
|---|---|---|
| S(2,2) | rank一non-square binary-quadratic tensor | rank二pencilのbasepoint-free型、common-factor型 |
| S(1,3) | pure-cube high block及びtransverse low blockのjoin型 | high binary-cubicの三相異根型、二重根+単根型 |

surface上のcentreは外点条件で除かれる。
表のACM帰属は一次資料Theorem6.2(b)(ii),(iii)、印刷17頁の全centre criterionを使う。
残る四orbitでは、全150二次係数を課してEuler商を取った十一parameter sourceを保持する。
十次元のquadric-times-constant-skew部分だけに制限していない。
conormal abelian必要条件の後の全二又は三parameter源で、actual Dは固定quadric因子を持つ。
全四model、全350 Jacobi及び126 actualD、correct projective bracket、extra sectionと全centre
網羅性は `PROJECTED_QUARTIC_SCROLL_SEPARATE_REVIEW.md` で別作者PASSである。
よって全四orbitの既約quintic枝は空である。

## 5. quartic coneからの全centre

Y=S(0,4)はrational normal quartic C4⊂P4のconeである。
centreをvertex shiftでbase P4へ移す。
base centreがclosed Sec(C4)内なら、一次資料Theorem6.2(a)、印刷17頁からimageはACM、従ってCI22。
closed Sec外ならC4の射影はsmooth rational quartic C⊂P3で、Sはそのconeである。

`PROJECTED_QUARTIC_CONE_RANK_ZERO_EXCLUSION.md` は全binary-quartic centre parameterを
一括して扱う。全Euler sourceを `Π=Q A+∂c∧V` と書き、nonCM extra Vを保持する。
その全sourceで `Π∧Π=Q*(Q*A∧A+2*A∧∂c∧V)` なのでactualDにはQ因子がある。
よって既約被約quintic枝は空である。
この弱い結論はJacobi、conormal条件、Qのrank、centre stabilizer又は全rank四不存在を
必要としない。cone noteの追加の強いproducer観察は本有限還元には使用しない。

## 6. Gへの論理と残る補集合

上記は全integral非退化degree四Sを覆い、CI22以外の枝を空、CI22枝をdimension≤24とする。
有限unionとPGL5-saturationのdimension boundにより、この支持を一般に持つ未知の
dimension≥25 Poisson componentは存在しない。
全CI22のknown-family帰属又は全source成分分類をこの弱いG排除から推論しない。

この結論はπ|S=0というrank零restrictionについてである。
rank二restriction、plane、及びdegree五・六の非normal支持は別の補集合として残る。
