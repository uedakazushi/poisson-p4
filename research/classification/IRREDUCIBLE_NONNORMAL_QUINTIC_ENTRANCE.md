# 既約非正規五次退化因子の曲面・導手・syzygyによる入口

2026-11-07。主義務 E の次数分割5、および G の未知優越成分を対象とする。
正規な実退化因子の全点分類は既存の `NORMAL_CLASSIFICATION` により閉じており、
ここでは再開しない。以下は**必要条件と条件付き認識**であって、
既約非正規五次をすべて T40 とする定理ではない。成分数は変更しない。

有限計算契約は `jobs/irreducible_nonnormal_quintic_entrance.json`、
実行器は `code/classification_nonnormal_quintic_entrance.py`、証明書は
`reports/CLASSIFICATION_IRREDUCIBLE_NONNORMAL_QUINTIC_ENTRANCE_20261107.json` である。
今回の幾何学的証明、一次資料からの入力、Macaulay2による確認を区別して記す。

## 1. 全既約非正規五次に対する有限数値入口

複素数体上で、\(D\subset\mathbb P^4\) を既約・被約・非正規な五次超曲面とする。
\(\nu:X\to D\) を正規化、\(L=\nu^*\mathcal O_D(1)\)、
\(M=\nu_*\mathcal O_X/\mathcal O_D\) と置く。
\(M\) の二次元既約支持を \(S_i\)、その一般点での長さを
\(\delta_i>0\) と書く。

**命題1（今回証明）。** 一般の平面切断の正規化の種数を \(g\) とすると

\[
 \sum_i\delta_i\deg S_i=6-g,\qquad 0\le g\le5. \tag{1}
\]

従って各 \(\deg S_i\le6\) であり、可能な**重み付き次数分割**は有限である。
これは曲面の同型類や埋込みの有限性を意味しない。
正規化のdivisorial conductorを \(C_X\ge0\) とすれば

\[
 K_X+C_X\sim0,\quad L^3=5,\quad
 L^2C_X=12-2g=2\sum_i\delta_i\deg S_i. \tag{2}
\]

**証明。** \(D\) はCohen--Macaulay三次元、\(\nu_*\mathcal O_X\) は \(S_2\) である。
深さ補題を各点で適用すると、\(M\) は支持上で局所Cohen--Macaulayな純二次元層となる。
例えば閉点では深さは少なくとも2、支持次元は2であり、
支持の曲線一般点・曲面一般点でも同じ局所次元に対応する深さを得る。
特に孤立した下位次元の寄与を隠していない。

一般の二超平面は \(M\) 上の正則列になる。平面五次曲線
\(C=D\cap H_1\cap H_2\) は既約であり、
\(\Gamma=X\cap\nu^*H_1\cap\nu^*H_2\) は \(\operatorname{Sing}X\) を避ける滑らかな曲線となる。
\(\Gamma\to C\) は有限双有理なので正規化である。切断後の完全列より

\[
 0\to\mathcal O_C\to\nu_*\mathcal O_\Gamma
       \to M|_C\to0
\]

を得る。右端の長さは \(\sum\delta_i\deg S_i\)、
平面五次の算術種数は6なので \(1\) が従う。

\(D\) はGorensteinで \(\omega_D\simeq\mathcal O_D\) である。
有限射の双対性より

\[
 \omega_X\simeq
 \mathcal Hom_D(\nu_*\mathcal O_X,\mathcal O_D)
       =\mathfrak c\mathcal O_X
       =\mathcal O_X(-C_X)
\]

をrank-one reflexive層として得る。これは \(K_X+C_X\sim0\) を意味する。
一般の \(\Gamma\) は \(X\) の滑らかな部分にあるので随伴公式から

\[
 2g-2=(K_X+2L)L^2=-L^2C_X+10
\]

となり \(2\) を得る。□

数値は \(g=0,1,2,3,4,5\) に対して、それぞれ

| \(g\) | \(\sum\delta_i\deg S_i\) | \(L^2C_X\) |
|---|---|---|
|0|6|12|
|1|5|10|
|2|4|8|
|3|3|6|
|4|2|4|
|5|1|2|

である。種数0または1への制限は、五次という次数だけからは出ない。

## 2. Poisson条件と、残すべきrank二つ

以下では \(D\) が大域Poisson bivector \(\pi\) の**実際の**退化因子であるとする。
Hamiltonian derivationは退化因子を保存し、従ってその正規性が破れる支持も保存する。
標数0では導分は正規化へ延長する。局所的に導手の元 \(a\) と正規化の元 \(b\) に対して

\[
 \partial(a)b=\partial(ab)-a\partial(b)
\]

であるから、導手も保存される。よって \(S_i\) はPoisson不変であり、
その一般点で \(\pi\) のrankは **0または2** である。
非正規性だけからrank0を課すことも、rank2を課すこともできない。

rank0なら \(\pi\in H^0(\mathbb P^4,E\otimes I_{S_i})\)、
ここで

\[
 E=\Lambda^2T_{\mathbb P^4}\simeq\Omega^2_{\mathbb P^4}(5).
\]

この線形条件はJacobi条件より弱く、次節の排除にはこれだけで十分である。

## 3. rank0曲面のACM排除と、非ACMの必須syzygy

**命題2（今回証明）。** rank0となる既約曲面 \(S\subset\mathbb P^4\) が
arithmetically Cohen--Macaulayなら、\(\deg S\le4\) である。

**証明。** Bottの公式、または外積Euler列から

\[
 H^1(E(-b))=H^3(E(-b))=0\quad(\text{全整数 }b),
\]

\[
 H^0(E(-d))=0\quad(d\ge3),\qquad
 H^2(E(-b))=
 \begin{cases}\mathbb C&b=5,\\0&b\ne5\end{cases} \tag{3}
\]

を得る。ACMなcodimension二の曲面idealはHilbert--Burch分解

\[
 0\to F_1\to F_0\to I_S\to0
\]

を持つ。\(3\)により \(H^0(E\otimes F_0)\to H^0(E\otimes I_S)\) は全射であり、
次数1または2のideal generatorだけが \(\pi\) に寄与し得る。

\(\deg S\ge5\) と仮定する。\(S\) が超平面に入るなら、
\(I_S=(\ell,f_d)\)、\(d=\deg S\ge5\) なので
\(\pi=\ell\sigma\) となる。非退化なら一次式はなく、独立な二つの二次式もない。
実際、\(I_S\) に入る二次式が可約なら素ideal性から一次式が入り矛盾する。
よって独立な二つの二次式は互いに素であり、Bézoutにより
\(\deg S\le4\) となってしまう。従って二次generatorは高々一つで、
\(\pi=f_2\sigma\)、または \(\pi=0\) となる。

どちらの非零の場合も実Pfaffianには \(\ell^2\) または \(f_2^2\) が因子として入り、

\[
 D=2V(\ell)+G_3\quad\text{または}\quad D=2V(f_2)+G_1
\]

となる。既約・被約な五次 \(D\) に反する。□

非ACMの場合も、自由分解のどの項がこの排除を回避し得るかを限定できる。

**命題3（今回証明）。** 飽和曲面idealの最小斉次自由分解を

\[
 0\to F_3\to F_2\to F_1\to F_0\to I_S\to0
\]

と書き、\(F_2\) に現れる \(\mathcal O(-5)\) の数を
\(\beta_{2,5}(I_S)\) とする。このとき

\[
 \dim\operatorname{coker}\big[
 H^0(EF_0)\to H^0(EI_S)\big]\le\beta_{2,5}(I_S). \tag{4}
\]

特に、rank0、\(\deg S\ge5\)、実 \(D\) が既約の場合は

\[
 \boxed{\beta_{2,5}(I_S)\ge1}. \tag{5}
\]

**証明。** \(K_1=\ker(F_0\to I_S)\)、
\(K_2=\ker(F_1\to K_1)\) とする。
三つの短完全列を \(E\) とtensorし、\(3\)を使うと、\(4\)左辺は
\(H^1(EK_1)\) に入り、さらに \(H^2(EK_2)\) に入る。
最後に \(H^2(EF_2)\to H^2(EK_2)\) は全射である。
従って左辺の次元は \(\dim H^2(EF_2)=\beta_{2,5}\) 以下となる。
\(\deg S\ge5\) では命題2の後半と同じBézout議論により
低次generatorの像は単一の一次式または二次式を因子に持つ。
従って既約 \(D\) の \(\pi\) はその像だけからは作れず、\(5\)を得る。□

これは非ACM曲面を排除する定理ではない。むしろT40の零曲面は、
まさに次数5の第二syzygyを使ってACMの場合の排除を回避する。

## 4. 一次資料におけるT40の正確な幾何

一次資料は Markarian--Polishchuk,
[*Compatible Poisson brackets associated with elliptic curves in G\(2,5\)*,
arXiv:2310.18759v2](https://arxiv.org/html/2310.18759v2) である。
今回§3.2のLemmas 3.2.1, 3.2.3, 3.2.4、Corollary 3.2.5と
§3.4のLemma 3.4.1\(i\)を直接読んだ。取得HTMLのhashと読解範囲は
`CLASSIFICATION_IRREDUCIBLE_NONNORMAL_QUINTIC_PRIMARY_20261107.json` に記録した。
以下の正規化・導手の結論は、引用した論文に存在しない全点分類を帰属させず、
記載された構造から導く。

楕円曲線 \(E_0\) とrank2・次数5のstable bundle \(V\) を取り、
\(\mathbb P^4=\mathbb P H^0(V)^*\) とする。\(q_{5,2}\) の退化五次を \(D_{E_0}\) と書く。
一次資料は次を証明している。

* 零locusの被約支持は、滑らかに埋め込まれた楕円scroll
  \(S_{E_0}=\mathbb P_{E_0}(V^\vee)\) である。
* rank二以下のlocus \(D_{E_0}\) は、\(A=\operatorname{Pic}^2(E_0)\) で添字付けられた
  平面 \(\mathbb P W_L\) の和である。それぞれの平面と \(S_{E_0}\) の交わりは
  \(E_0\) の次数3の埋込みである。
* rank3bundle \(W\to A\) による射
  \(f:X_0=\mathbb P_A(W)\to D_{E_0}\) は、\(S_{E_0}\) の外では同型である。
* \(S_{E_0}\) の点に対応するextensionのmiddle bundleは
  \(\mathcal O(p)\oplus V'\)、\(V'\) はrank2・次数4のsemistable bundleで
  \(\dim\operatorname{End}(V')=2\) である。
  \(V'=L_1\oplus L_2\)（非同型の次数2線bundle）、または非分解型である。

**命題4（上記一次資料から今回導出）。** \(f\) は有限であり、従って \(D_{E_0}\) の正規化である。
その導手divisorは滑らかな既約曲面

\[
 B\simeq E_0\times\operatorname{Pic}^2(E_0)\subset X_0
\]

であり、各平面fiberでは次数3曲線となる。
\(B\to S_{E_0}\) は一般に二対一である。
導手idealは \(X_0\) 上では \(\mathcal O_{X_0}(-B)\)、
\(D_{E_0}\) 上では被約ideal \(I_{S_{E_0}/D_{E_0}}\) となる。

**証明。** \(S_{E_0}\) の点上の \(f\) のfiberは、
そのextensionがsplitする次数2線subbundleの選択である。
次数2線bundleから \(\mathcal O(p)\) への射は零なので、それは \(V'\) に入る。
\(V'=L_1\oplus L_2\) なら選択はちょうど二つ。
非分解型なら \(V'\) は次数2線bundle \(M\) の非自明な自己extensionであり、
次数2線subbundleは一つである。従って全fiberは有限集合で、
properな \(f\) はquasi-finite、従って有限となる。
\(X_0\) は滑らかで \(f\) は双有理なので正規化である。

各平面内の次数3の \(E_0\) は相対的に同じ曲線の族を成し、
その全体 \(B=E_0\times A\) は \(f^{-1}(S_{E_0})\) の被約支持である。
一般点での二つの選択を交換するinvolutionは

\[
 (p,L)\longmapsto
 (p,\det(V)\otimes\mathcal O(-p)\otimes L^{-1}). \tag{6}
\]

固定条件 \(L^2=\det(V)(-p)\) は各 \(p\) に対して四点を持つ。
ここから、特殊点の全解析的標準形までは主張しない。

Lemma 3.4.1の一般点での線形化は

\[
 \pi_{\rm lin}=2e_1\partial_{h_1}\wedge\partial_{e_1}
 -e_1\partial_{h_2}\wedge\partial_{e_1}
 +2e_2\partial_{h_2}\wedge\partial_{e_2}
 -e_2\partial_{h_1}\wedge\partial_{e_2}.
\]

その全Jacobiが零、Pfaffianが \( -3e_1e_2\) であることを今回も展開した。
従って \(S_{E_0}\) の一般点を横切る曲線特異点はordinary nodeで、\(\delta=1\) である。
\(X_0\) は滑らかなのでその導手は可逆層であり、
支持は既約な \(B\)、一般点での位数は1である。従って導手divisorはちょうど \(B\) となる。
さらに局所的に \(\mathcal O_D/\mathfrak c\hookrightarrow
\nu_*\mathcal O_{X_0}/\mathfrak c\mathcal O_{X_0}\) で右辺は被約である。
よって \(D\) 上の導手も被約で、支持が \(S_{E_0}\) なのでそのidealに等しい。□

tautological交点数から \(\deg S_{E_0}=5\) である。
命題1により \(g=1\)、\(L^2B=10\) となる。
\(S_{E_0}\) 外では正規化が滑らかな同型であり、
\(S_{E_0}\) 上では \(\pi=0\) なので五次式の一次微分も零となる。
従って \(\operatorname{Sing}(D_{E_0})_{\rm red}=S_{E_0}\) である。

## 5. 正規化上のresidueから言えることと言えないこと

\(\omega=\pi^{-1}\) は \(\mathbb P^4\setminus D\) の閉二形式である。
\(D_{\rm reg}\) では対数極を持ち、非零の閉一形式

\[
 \alpha=\operatorname{Res}_{D_{\rm reg}}\omega
\]

を与える。これを \(X\) 上の有理閉一形式と見ると、divisorialな極の支持は導手に含まれる。
ただし、**正規化後も極が高々一位であるとは限らない**。

具体的な局所反例は

\[
 \pi=\partial_s\wedge\partial_t+
       (v^2-u^3)\partial_u\wedge\partial_v.
\]

これは全Jacobiを満たし、\(D=(v^2-u^3=0)\) は既約なcuspidal因子である。
\(u=a^2,v=a^3\) による正規化では

\[
 \alpha=\frac{du}{2v}=\frac{da}{a^2}=-d(1/a)
\]

となり二位の極を持つ（residueの全体符号は規約に依存する）。
この例は**局所Poisson germ**であり、
大域的な \(\mathbb P^4\) の反例を構成したとは主張しない。

**命題5（条件付き、今回証明）。** \(\alpha\) のdivisorialな極が高々一位であり、
導手の被約支持 \(B\subset X\) が既約だと仮定する。
このとき \(\alpha\) は非零の閉reflexive一形式となり、

\[
 g\ge1,\qquad \sum_i\delta_i\deg S_i\le5. \tag{7}
\]

**証明。** \(X\) と \(B\) の一般滑らかな点で、閉性と一位の極により
\(\alpha\) はlogarithmicとなり、そのresidueは定数 \(\lambda\) である。
一般の \(\Gamma=L_1\cap L_2\) を取ると、\(\alpha|_\Gamma\) の極は
\(\Gamma\cap B\) に限り、その全residue和は
\(\lambda L^2B=0\) である。
\(L\) はample、\(B\) は非零有効なので \(\lambda=0\) となる。
従って全codimension一で極が消え、正規性から
\(\alpha\in H^0(X,\Omega_X^{[1]})\) を得る。
一般の \(\Gamma\) は非零形式 \(\alpha\) の核に接しないので、
非零の正則一形式 \(\alpha|_\Gamma\) を持つ。従って \(g\ge1\) である。□

\(X\) が特異な場合、reflexive一形式がresolutionへ正則に延長するとは
この議論だけからは結論できない。従って一般に \(q(\widetilde X)>0\) とは書かない。
\(X\) が滑らか、またはこの形式の正則延長を別に証明した場合には、
非零正則一形式により \(q(\widetilde X)>0\) となる。
一般ordinary-double型では一位極の仮定は局所計算で検査できるが、
cuspidal型を同時に覆う仮定ではない。

ordinary-double型の連結な導手二重被覆では、枝交換は二つのresidueを符号反転させる。
これによりlogarithmic二形式の \(du/u\wedge dv/v\) の係数が零になっても、
**Poisson構造の不存在は従わない**。横断lognormal方向と接方向の交差項で
非退化性を保つことができ、その場合は二重曲面上のrank0へ送られる。
実際のT40はこの可能性を保持しなければならない例である。

一方、導手が \(B_1+B_2\) にsplitし、\(X\) が滑らかで、
この導手がSNC、かつ対数形式の仮定が成り立つ場合、
residueが \(\lambda,-\lambda\) ならGysin関係は

\[
 \lambda([B_1]-[B_2])=0\quad\text{in }H^2(X,\mathbb C). \tag{8}
\]

\(\lambda\ne0\) で得るのはまずcohomology classの一致である。
Picard零次成分やtorsionを排除せずに線形同値へ読み替えない。
追加条件で \(B_1\sim B_2\) を証明したときだけ、
\(-K_X\sim2B_1\) と対応するpencilの議論へ進める。

## 6. generic T40零scrollの正確なsyzygyと条件付き認識

既存 `code/elliptic_weight1_dominance.py` と同じcyclic \(q_{5,2}\) の斉次五変数表示を使い、
**標本点ではなく \(\mathbb Q(\ell)\)** 上で計算した。
全十個のraw Jacobiを零と確認した後、
Euler場 \(\mathcal E\) による十個の三次式 \(\mathcal E\wedge\pi\) を作った。
射影零locusはこのidealの飽和であり、得たideal \(I_S\) について

\[
 \operatorname{codim}S=2,\quad\deg S=5,\quad
 0\to\mathcal O(-5)\to\mathcal O(-4)^5
 \to\mathcal O(-3)^5\to I_S\to0 \tag{9}
\]

をMacaulay2 1.24.11で確認した。
全differentialの積、Buchsbaum--Eisenbudの必要grade
\(\operatorname{codim}I_4(d_2)\ge2\)、
\(\operatorname{codim}I_1(d_3)\ge3\) も検査した。
さらに \(I_S+I_2(\operatorname{Jac} I_S)\) の飽和が単位idealとなるので、
generic零曲面の射影滑らかさも確認された。
全M2入力と全出力行列は証明書に埋め込んで保存した。
実行約3.20秒、Python controllerの最大RSSは12,672 KiB、
子プロセスを含むアドレス空間上限は1,536 MiBである。
このRSS値をM2全体の使用量としては報告しない。
独立した第二CAS実装による再実行は行っていない。

**系6（今回証明、対象scrollを特定する条件付き）。** \(9\)を持ち、
既知の \(q_{5,2}\) bivectorが消える同じ埋込み曲面 \(S\) に対して

\[
 H^0(\mathbb P^4,E\otimes I_S)=\mathbb C\pi_{5,2}. \tag{10}
\]

実際、\(9\)の \(F_0\) は三次generatorだけを持つので \(H^0(EF_0)=0\) であり、
命題3から左辺次元は高々1である。既知非零 \(q_{5,2}\) が存在するので次元は1となる。
従って、ある未知bivectorのrank0曲面を、\(9\)を持つ**この既知scroll**と
射影同値に特定できれば、そのbivectorはT40に入る。

これは「同じBetti数の全曲面が既知scrollである」とも、
「非ACMな次数5または6曲面はすべてそのscrollである」とも主張しない。
今回の有理関数体上の計算はgeneric parameterの証明であり、
特殊な楕円parameterでの全点resolutionについては別の幾何学的証明が必要である。
元の有限計算契約にもこの補集合を残した。

## 7. 今回の還元後に残る正確な課題

既約非正規五次の全支持を閉じるため、少なくとも次を保持する。

1. 非正規曲面の一般rank2の全配置。重み付き次数は6以下だが、
   rank0用のideal-generator排除は適用しない。
2. rank0で次数4以下のACM曲面。
3. rank0で次数6以下の非ACM曲面。次数5または6なら
   \(\beta_{2,5}\ge1\) が必要であるが、それだけではT40と認識できない。
4. 正規化上のresidueが高位極を持つ型、および導手の被約支持が可約な型。
   命題5の仮定を満たす補集合だけを暗黙に捨てない。
5. \(\beta_{2,5}\) による残る曲面の分類、またはPoisson線形化と導手の追加条件から
   既知楕円scrollへの認識を導く命題。

未知G優越成分の射影次元が25、26、27に限られることは
`DOMINANT_COMPONENT_GEOMETRY.md` の既存結果である。
上の曲面次数制限だけから、そのstratumの全PGL飽和次元やprimitive一次底への
非優越性は結論できない。曲面と許されるbivectorのincidence次元を別に評価する必要がある。
今回の具体的な不足は、\(5\)を満たす小次数非ACM曲面と
rank2曲面を、実Poisson条件・導手条件の下で尽くす認識定理である。
