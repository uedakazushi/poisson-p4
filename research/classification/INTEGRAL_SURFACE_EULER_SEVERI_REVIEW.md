# integral曲面のSeveri・Euler補題：特異曲面の射程を含む独立読解

2026-11-07。次の二稿を独立に監査し、記載された仮定と結論を受理する。

* INTEGRAL_SURFACE_LINEAR_NORMALITY.md
* RANK_ZERO_REGULAR_LINEARNORMAL_SURFACES.md

新しい係数計算・source probe・degree3計算は行っていない。
一次資料の読解場所とinput hashは
CLASSIFICATION_INTEGRAL_SURFACE_EULER_SEVERI_REVIEW_20261107.json に記録した。

## 1. 一次資料はsmooth版だけの定理ではない

[Chiantini, Lectures on the structure of projective embeddings](https://emis.de/ft/47035)
を既設web readerで直接開いた。
printed p.335の全体設定は複素数体上の既約・被約・非退化projective varietyである。
printed pp.369–370、Theorem14は1-defective surfaceをconeまたはVeroneseへ分類し、
当該sectionにも追加のsmoothness仮定はない。

独立に [Niu–Park, arXiv:2203.04601v1](https://arxiv.org/pdf/2203.04601v1) の
PDFを取得して読み、printed pp.15–16、Proposition5.3の証明冒頭を確認した。
そこでは任意の非退化surface in P4について、同型射影されたVeronese以外は
linearly normalであるという補題を先に証明し、特異曲面版Severiとconeのvertex tangentを使う。
同propositionのdimension三以上というstatementを、曲面へ直接誤適用していない。
Dale原論文そのものの独立精読は今回行っていない。

## 2. Hilb²を使うclosed secantへの橋は全長さ二schemeで成立する

Sを任意のintegral・非退化surface in P4、\(H=\mathcal O_S(1)\) とする。
\(h^0(H)=r+1>5\) なら完全系は閉埋込み \(j:S\hookrightarrow\mathbb P^r\) を与える。
元の埋込みは五次元subsystem Wによる、中心
\(\Lambda=\mathbb P^{r-5}\) からの射影である。

Wが閉埋込みを与えるので、**特異点に支持されたものも含む**全長さ二scheme Zに対して
\(W\to H^0(Z,H|_Z)\) は全射である。
従って \(\Lambda\) は各span line \(\langle j(Z)\rangle\) を避ける。
Hilb²(S)はproper、普遍長さ二族のtautological bundleはrank2で、
完全系のevaluationが全点で全射なのでGrassmannへのspan mapは全域で定義される。
その像上の普遍lineの全体もproperであり、P^rへの像は閉集合となる。
その集合は相異なる二点の全secantを含み、従ってclosed secant varietyを含む。

よって \(\Lambda\cap\operatorname{Sec}(j(S))=\varnothing\)。
Hilb²の滑らかさや既約性を仮定する必要はない。
projective dimension theoremから \(\dim\operatorname{Sec}(j(S))\le4\) となり、
r≥5での期待次元5より小さい。上記の特異Severi定理を適用できる。

## 3. cone例外はsmoothnessではなくembedding dimensionで除ける

j(S)がconeなら、そのvertexの点pと任意の点qを結ぶlineがj(S)に入る。
j(S)は非退化なので、それらの方向は全 \(T_p\mathbb P^r\) を張る。
従ってpのintrinsic embedding dimensionはr≥5である。
一方、同じ抽象scheme SはP4へ閉埋込みされているため、その局所embedding dimensionは4以下。
これは矛盾であり、coneを除ける。

残る完全埋込みは \(v_2(\mathbb P^2)\subset\mathbb P^5\)。
従ってlinearly normalでないintegral非退化surface in P4は、
完全Veroneseからの同型射影であり、抽象的には滑らかなP2、HはO(2)である。
singular coneの正規化に余分なsectionがあることを、
coneそのもののlinear nonnormalityと読み替えていない。

## 4. restricted Euler snakeはnormality・ACM性を使わない

\(V=\mathbb C^5\)、\(\mathcal E=\Lambda^2T_{\mathbb P^4}\) と書く。
外積Euler列は

\[
 0\to T_{\mathbb P^4}\to\Lambda^2V\otimes\mathcal O(2)
                \to\mathcal E\to0.
\]

全てvector bundleなのでSへのrestrictionもexactである。
\(H^1(T_{\mathbb P^4})=0\) により任意のglobal bivectorは全150二次係数のliftを持つ。
\(\pi|_S=0\) なら、そのliftのrestrictionは一意の
\(\tau\in H^0(S,T_{\mathbb P^4}|_S)\) の像である。
liftをEuler-vertical項で変えた場合、tauはambient線形場のrestrictionだけ変わる。
従って

\[
 0\to I_S(2)\otimes\Lambda^2V
  \to H^0(\mathcal E I_S)
  \to\operatorname{coker}\!\left[H^0(T_{\mathbb P^4})
                    \to H^0(T_{\mathbb P^4}|_S)\right]
\]

がexactとなる。右端への全射は主張しておらず、上界には不要である。

左端のinjectivityを確認する。
I_S係数のliftがEuler-verticalなら、その線形projective場はS上で零である。
任意のmatrixのprojective fixed locusは有限個のeigenspaceの和であり、
integral Sは一つに含まれる。SがP4をspanするのでmatrixはscalarで、vertical bivectorも零。
この議論はmatrixの対角化可能性を仮定しない。
同時に、ambient線形場のrestrictionのinjectivityも証明される。

restricted通常Euler列

\[
 0\to\mathcal O_S\to\mathcal O_S(1)^5\to T_{\mathbb P^4}|_S\to0
\]

から、\(r_2=h^0(I_S(2))\)、\(l=h^0(\mathcal O_S(1))-5\)、\(q=h^1(\mathcal O_S)\) として

\[
 h^0(\mathcal E I_S)\le10r_2+5l+q.
\]

q=l=0ならambient線形場のrestrictionは同型となり、完全sourceは正確に

\[
 H^0(\mathcal E I_S)=I_S(2)\otimes\Lambda^2V.
\]

追加の \(H^1(\mathcal O_S(1))\) 消滅、normality、smoothness、ACM性は不要である。

## 5. degree五以上・q=0の排除の範囲

Severi例外のprojected Veroneseはdegree四である。
従ってintegral・非退化・degree五以上のSはlinearly normalであり、
q=0の場合は上の完全quadratic sourceとなる。
独立な二quadricがSを含めば、非退化性と素ideal性からそれらは互いに素で、
degree四のcomplete intersectionにSが含まれることになる。従って \(r_2\le1\)。

\(r_2=0\) なら全bivectorが零。
\(r_2=1\) ならliftは \(Q_2 A\)、Aはconstant skew matrixであり、実五次退化式は

\[
 D=Q_2^2\,\iota_E(A^2/2)
\]

すなわちQ²と一次式の積、または恒等零となる。
Euler-vertical補正は \(E\wedge\Pi^2\) を変えないので、この計算は実退化因子に対するもの。
従って既約・被約五次を持つrank4 Poisson構造のrank0曲面にはなれない。
Jacobiを課す前の排除である。

normal Sでresolutionのqが零ならLerayの単射からq(S)=0を得るが、
非normal Sの正規化がrationalであるだけではq(S)=0と仮定しない。
q>0、degree四以下の他の軌道、rank2曲面はこの補題では残る。
既約非正規五次の全分類、全Gの完成、新成分数の変更は結論しない。
