# smooth quartic scroll の四 projection orbit：独立監査

2026-11-07、主義務 G。**PASS**。
`PROJECTED_QUARTIC_SCROLL_ISOLATED_MODELS.md` の四非ACMモデルについて、
全 projective rank零・既約被約 quintic 枝は空である。
残る smooth S(2,2),S(1,3) の中心は実際に ACM→CI22 へ接続され、捨てた中心はない。
S(0,4) cone の一 parameter と衝突型は別補集合として保持する。
全 degree四曲面又は G 全体の完成をこの監査から結論しない。

## 1. 中心の全 orbit と ACM への接続

S22 の六座標は `(su²,suv,sv²,tu²,tuv,tv²)`。
中心を C²⊗Sym²C² の rank で分ける。
rank1 で binary quadratic が平方なら中心は曲面上で、投影の定義域外。
rank1 nonsquare は Segre P1×P2 上の残りであり、下記一次定理により ACM image へ入る。
rank2 は二次元 binary-quadratic pencil の row space を定める。
その annihilator は P(Sym²C²*) の点で、base PGL2 の orbit は rank2 と rank1 の二つ。
従って pencil は `(u²,v²)` 又は `(u²,uv)` に直せる。
fiber PGL2 が全 row basis を取り替えるため、これは実際の二つの whole centre orbit。
代表中心と quotient が S22 node/cusp の二モデルである。

S13 の六座標は `(Us,Ut,Vs³,Vs²t,Vst²,Vt³)`。
中心の cubic block C=(c0,c1,c2,c3) が零なら directrix 上で除外される。
bundle automorphism `U→U+V(αs²+βst+γt²)` は実際の ambient linear map を誘導し、
linear block を

\[
 (a_0,a_1)\mapsto
 (a_0+\alpha c_0+\beta c_1+\gamma c_2,
  a_1+\alpha c_1+\beta c_2+\gamma c_3)
\]

と変える。Hankel 行列の rank は C が pure cubic power の場合だけ1、その他は2。
Sym³ のこの Veronese 座標を通常の polynomial coefficient に直す際は binomial 係数を含む同一視を取る。
rank2 なら linear block を全て零に出来、残る cubic は三相異根又は二重根+単根の二型。
base PGL2 はそれぞれに推移的なので、C=(1,0,0,1) と C=(0,1,0,0) の二モデルで全てを覆う。
pure C なら linear block の一成分を shear で消す。
残りも零なら曲面上、非零なら `Join(S(1),S13)\S13` に入る一つの orbit であり ACM image となる。

この ACM 判定は Brodmann–Park,
*On varieties of almost minimal degree I: Secant loci of rational normal scrolls*,
[arXiv:0808.0090v3](https://arxiv.org/pdf/0808.0090v3) の Theorem6.2(b)(ii),(iii)、印刷17頁と
Remark6.4(v)、印刷19頁を独立に直接確認した。
codim≥2、代数閉体、中心が元曲面外という同定理の全仮定を満たす。
最初の UZH URL は timeout だったが、arXiv の同一次資料を開いて確認した。

中心外からの投影は finite：正次元 fiber なら中心を通る直線全体が元曲面に入り矛盾する。
image は非退化で、次数の積 formula と degree4 から birational・degree4 である。
次数2以下の非退化曲面は P4 に存在しないため、次数を落とす多重 cover はない。
ACM image の三つの linear parameters による Artin quotient は embedding dimension2、length4、
Hilbert function(1,2,1)。従って I2 は二次元。
素性・非退化性から二 quadric は共通一次因子を持たず、degree4 の完全交叉を与える。
degree が image と一致し、完全交叉は unmixed なので image 自身が scheme-theoretic CI22 である。

この枝は、別に受理済みの regular-pencil V2 と
`ALL_DEGENERATE_CI22_SEPARATE_REVIEW.md` の全退化 pencil 証明を合わせた
**全 integral 非退化 CI22 の rank零 G 次元≤24**へ接続する。
この既済依存を四非ACMモデルの不存在証明と混同しない。

中心が normalization の実 automorphism で移れば、六空間の中心直線商の五次元線形同型を誘導する。
従って各モデルの source 計算は一標本ではなく、その whole PGL5 surface orbit に適用する。
逆に全中心は上の rank/根型の分岐を通るので、四モデル又は ACM image 以外の未記載中心はない。

## 2. 全150から Euler 商十一まで

独立実装は `code/classification_projected_quartic_scroll_review.py`。
producer、SymPy、CAS を import せず、別作者の Fraction 疎多項式 helper で再構成した。
全150係数の quadratic bivector Π を保持する。
外積 Euler と H1(TP4)=0 により、これが全 global projective section を表し、kernel は
Euler∧linear field の二十四次元である。

normalization parametrization q に対し、projective zero は正確に `q∧Π(q)=0`。
raw Π(q)=0 という強い別条件を使わない。
S22 の全280行、S13 の weighted 全300行の150列行列を、全十 triple と全 monomial から作った。
公開 matrix と全係数一致し、rank115 minor と完全 kernel35を検証した。
Euler rank24、全行列への積零も再生した。

canonical projector `Π→Π−E∧divΠ/5` により、完全150×11 basis の発散が零であることと、
全 kernel35 の canonical 像がこの十一に入ることを確認した。
Euler24 と十一 basis の直和 rank は35であり、全 quotient が得られる。
独立に quadric 制限行列の rank14を確認して I2=CQ を再生し、
canonical な `Q×constant skew` 十クラスの rank10と、十一 basis の中での包含も確認した。
**nonCM image の extra 一クラスを含む全十一を残している。**
ACM と仮定して Q×十クラスだけに縮めていない。

## 3. generic conormal と正しい projective 補正

四組 F2,F3 の pullback は零で、Jacobian の2×2 minor の少なくとも一つが pullback 後に非零。
各 parametrization は元の scroll の dense chart を写す dominant map であり、
generic に F2,F3 は S の normal basis をなす。
S22 node は x2x3、S22 cusp は x2x4、S13 は x3x4 を反転した base ratio の式からも確認出来る。
global V(F2,F3) の residual component を S と同一視していない。

既約 quintic と rank零仮定から、非plane secant/Frobenius 補題により generic conormal は abelian。
x0≠0 上の関数 F2/x0²、F3/x0³ の bracket を展開すると、その numerator は厳密に

\[
 K=x_0\{F_2,F_3\}-3F_3\{F_2,x_0\}-2F_2\{x_0,F_3\}.
\]

実際の bracket は K/x0⁶。
raw Π が S 上に Euler-vertical 値を持つ可能性を残すので、この二つの補正項を落とさない。
abelianity は generic に K∈I² を必要とし、従って全 ambient first partial は S 上で零。
この必要方向だけを用い、nonCM な global ideal に ordinary-square/complete-intersection 論を誤用しない。
generic chart の反転は条件導出だけに使い、その後の polynomial pullback 零性は integral S 全体へ延びる。
chart の外の別 Poisson 源を捨てた計算ではない。

全十一 basis について K 自体と五つの全 first partial を作り直し、
S22 の全330×11、S13 の全360×11行列を公開データと一致させた。
rank minor と完全 kernel を検証した。

|モデル|全 raw kernel|Euler 商|conormal rank|全必要 source|
|---|---:|---:|---:|---:|
|S22 node|35|11|9|2|
|S22 cusp|35|11|9|2|
|S13 node|35|11|9|2|
|S13 cusp|35|11|8|3|

## 4. 全二・三 parameter 源の fixed Q 因子

完全 kernel を十一 basis に戻した全150×2 又は150×3 matrixを再生し、
全自由 parameter z0,z1,z2 を残した実 tensor を構成した。
全350 Jacobi polynomial、parameter 二次係数行列、その rank1,1,1,3、
全126 signed Euler-Pfaffian 係数をそれぞれ独立に一致させた。

その actual D を Q で独立な疎多項式長除算にかけ、全四型で余り零、
商が `CLASSIFICATION_PROJECTED_QUARTIC_SCROLL_FACTORS_20261107.json` の cubic と一致することを確認した。
商の疎係数は独立証明書にも保存する。

\[
 D_{\rm actual}=Q\cdot C_3(z_0,z_1,z_2;x_0,\ldots,x_4)
\]

は **Jacobi 前の全 parameter 源**に成り立つ。
従って全必要 Poisson 源の actual quintic は零又は fixed quadric 因子付きであり、既約被約 quintic ではない。
源の selected point、source parameter の主開部分、数値標本、GB 結果への推論を使わない。

S13 の z1²、cusp の追加 z2²についても、保存された元350 Jacobi 行の functional を直接積み戻した。
これらが零となる全被約支持上で D=0 になることも確認したが、G 排除には Q 因子だけで十分である。

## 5. 実行・保存・未処理補集合

有限契約 `reports/CLASSIFICATION_PROJECTED_QUARTIC_SCROLL_REVIEW_JOB_20261107.json` の gate PASS 後、
独立 checker の全再構成は1.9611秒、最大 RSS21,564KiB、全 PASS。
結果は `reports/CLASSIFICATION_PROJECTED_QUARTIC_SCROLL_REVIEW_20261107.json`。
全入力 hash、確定数学ノート snapshot、checker 全文、全 source rank minor、独立 cubic quotient を保持する。
初回 checker 内の helper の binary multiplication 呼出を訂正した後の最終実行であり、
producer・旧証拠 bytes は変更していない。

今回の結論は四非ACM whole orbit の不可約 quintic 枝の空性と、他の smooth-scroll 中心の ACM/CI22 接続。
S(0,4) の投影 family、一 parameter、及び衝突型をこの四モデルで置き換えない。
新成分、全 E 帰属、G 全体の完成、全 degree四曲面の排除は主張しない。
共有台帳・原稿・Git は編集していない。
