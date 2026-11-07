# 滑らかなrank0曲面のclean ordinary-double枝をT40へ認識する

2026-11-07。既約非正規五次の主分岐に対する認識定理である。
先のordinary-double認識で課した「正規化Xが滑らか」「導手Cが滑らかで既約」という
二仮定を用いず、実際の非正規曲面Sとそのblowupを使う。
verticalな例外交叉、特異なS、rank2、他の横断特異型は捨てない。

主たる局所pole議論と小次数曲面補題はroot提案を独立に読解し、
以下の証明として整理した。有限計算契約は
jobs/irreducible_nonnormal_clean_rank_zero.json、
実行器は code/classification_nonnormal_quintic_entrance.py の --clean-rank-zero、
証明書は CLASSIFICATION_IRREDUCIBLE_NONNORMAL_CLEAN_RANK_ZERO_20261107.json。
gateと恒等式算術は幾何学的仮定の成立自体を検証するものではない。

## 1. 正確な仮定と結論

\(\pi\) を複素 \(\mathbb P^4\) のrank4 Poisson構造、
\(D\) をその実際の既約・被約・非正規五次退化因子とする。
\(S\subset\operatorname{Sing}D\) を非正規locusの二次元既約成分とし、次を仮定する。

1. \(S\) は滑らか。
2. \(\pi|_S=0\)。一般点でrank0という仮定でも、正則性と被約性からこの条件を得る。
3. \(S\) の一般点で、Dは横断ordinary doubleである。
4. \(b:Y=\operatorname{Bl}_S\mathbb P^4\to\mathbb P^4\)、
   \(p:E=\mathbb P(N_{S/\mathbb P^4})\to S\)、strict transform \(\widetilde D\) とすると、
   \(\Gamma=\widetilde D|_E\) は**verticalな既約divisor成分を持たない**。
   すなわち全既約成分がSをdominantに写す。

仮定4をclean条件と呼ぶ。この条件は一般ordinary-double性だけから自動とはしない。
E上でvertical divisorとなるのは、Sの曲線上の全 \(\mathbb P^1\) fiberを含む場合である。
Gammaの分岐、特異点、水平成分同士の交叉を滑らかと仮定する必要はない。

**認識定理（今回証明）。**
この四仮定の下でSは次数5の滑らかな楕円scrollであり、
ある楕円曲線Bとrank2・次数5 stable bundle Fについて

\[
 (S,\mathcal O_S(1))\simeq
  (\mathbb P_B(F^\vee),\mathcal O_{\mathbb P(F^\vee)}(1)).
\]

\(\pi\) はこの同じ埋込みscrollのcanonical \(q_{5,2}\) のscalar倍である。
従ってT40に属す。他の非正規成分や正規化の全域を先に分類する必要はない。

## 2. rank0によりblowup上の逆二形式のE極は高々一位

\(\omega=\pi^{-1}\) は \(\mathbb P^4\setminus D\) 上の閉二形式である。
四次元のskew行列のinverseは、各係数が相補2×2係数をPfaffianで割った形となる。
従って、Sの局所idealを \(I_S\) とすると、
rank0条件によりinverseの全分子は \(I_S\) に入り、分母はDの局所方程式である。

一般ordinary-double性から \(\operatorname{ord}_S D=2\) である。
blowupのE方程式をtとすると、分子はtで割れ、分母は
\(t^2\widetilde f\) となる。微分形式のpullbackは負のt次数を加えないので、

\[
 b^*\omega\in
 \Omega_Y^2\bigl(*\widetilde D+E\bigr)
\]

で、Eに沿う極は高々一位である。
ここで正規化の滑らかさ、導手の成分数は使っていない。

閉性からEの一般点で \(b^*\omega\) はlogarithmicとなり、
有理閉一形式

\[
 \beta=\operatorname{Res}_E b^*\omega
\]

を得る。betaの候補となるdivisorial poleは正確にGammaの支持に含まれる。
実際、Eでresidueを取る前の係数の分母は \(\widetilde f\) のみであり、
Eへ制限するとその零divisorはGammaである。
分岐discriminantの因子を、別の分母として付け加えていない。

## 3. 全水平poleが消え、betaは非零正則形式となる

Gammaの各水平成分の一般点はSの一般ordinary-double点上にある。
そこでDの二枝を \(uv=0\)、Sを \(u=v=0\)、
接方向を \(z,w\) とする。必要ならこの局所計算だけétaleに二枝をsplitする。
log基底 \(du/u,dv/v,dz,dw\) によるomegaの行列を

\[
 W=\begin{pmatrix}
 0&\lambda&a&b\\
 -\lambda&0&c&d\\
 -a&-c&0&e\\
 -b&-d&-e&0
 \end{pmatrix},\qquad
 \operatorname{Pf}W=\lambda e-ad+bc
\]

と書く。実Dが被約ordinary doubleであることから
\(\operatorname{Pf}W\) はこの一般点でunitである。
実bivectorの接方向成分は
\(-\lambda(\operatorname{Pf}W)^{-1}\partial_z\wedge\partial_w\)。
よってrank0は \(\lambda|_S=0\) を強制する。

blowup chart \(u=t,v=ts\) では

\[
 \frac{du}{u}=\frac{dt}{t},\qquad
 \frac{dv}{v}=\frac{dt}{t}+\frac{ds}{s},
\]

従ってE上のresidueは

\[
 \beta=\lambda|_S\,\frac{ds}{s}
       +(a+c)|_S\,dz+(b+d)|_S\,dw. \tag{1}
\]

水平Gammaの一位極のresidueは零であり、betaはその一般点で正則である。
もう一つのblowup chartも同じで、全水平成分のpoleが消える。
clean条件により他のvertical divisor成分はない。
従ってbetaはEの全codimension一で正則である。
Eは滑らかなので \(\Omega_E^1\) はreflexiveであり、Hartogsにより

\[
 \beta\in H^0(E,\Omega_E^1)
\]

となる。Gammaの分岐点や特異点はcodimension二以下の例外集合としてこの延長に含まれる。
それらを無断で除外したのではない。

しかも \(\lambda|_S=0\) でのunit性は \(ad-bc\ne0\) を意味する。
もし \((a+c,b+d)=(0,0)\) なら \(ad-bc=0\) となるので、(1)は非零である。
従ってbetaは非零の大域正則一形式である。

projective bundleの相対cotangent列と
\(H^0(\mathbb P^1,\Omega_{\mathbb P^1}^1)=0\) から

\[
 H^0(E,\Omega_E^1)=p^*H^0(S,\Omega_S^1).
\]

よって

\[
 \boxed{q(S)=h^0(S,\Omega_S^1)\ge1}. \tag{2}
\]

式(1)の全log-skew行列変換、非零性のcoframe恒等式は有限算術でも確認した。
poleの支持、Hartogs、projective bundleからの降下は本文の幾何学的証明である。

## 4. 小次数の滑らかなirregular曲面に関する補題

**補題。** 滑らかな既約曲面 \(S\subset\mathbb P^4\) が
\(q(S)\ge1,\ d=\deg S\le6\) を満たすなら、次数5の楕円scrollである。

以下で証明する。Hはhyperplane divisor、
\(g_H=1+(H^2+K_SH)/2\) は一般hyperplane曲線の種数とする。
Sが \(\mathbb P^3\) に入ればhypersurfaceの完全列より \(q(S)=0\) であるから、
Sは非退化である。一般hyperplane曲線も \(\mathbb P^3\) の非退化曲線となる。
Castelnuovo上界は

\[
 g_H\le\left\lfloor\frac{(d-2)^2}{4}\right\rfloor.
\]

従って \(K_SH=2g_H-2-d<0\) であり、唯一の例外候補は
\(d=6,g_H=4,K_SH=0\) である。

もし \(\kappa(S)\ge0\) なら非零pluricanonical sectionがあり、
ampleなHとの交点から例外候補以外は不可能である。
例外候補ではその有効divisorのH次数が零なので、
あるmについて \(mK_S\sim0\) となる。
従って \(K_S^2=0,\ p_g(S)\le1,\ \chi(\mathcal O_S)\le1\)。
一方、滑らかな曲面の \(\mathbb P^4\) double-point formula

\[
 d^2-10d-5K_SH-2K_S^2+12\chi(\mathcal O_S)=0 \tag{3}
\]

はこの例外で \(\chi(\mathcal O_S)=2\) を強制し、矛盾する。

よって \(\kappa(S)=-\infty\)。
ここでは標準的な滑らかな射影曲面のruledness定理
「Kodaira次元 \(-\infty\) ならrational又はgeometrically ruled surfaceのblowup」
を使う。irregularなので、底の種数 \(g=q(S)\ge1\) のruled surfaceのblowupである。
未確認のpolarized三次元多様体分類やFujita定理は使わない。

## 5. ruled surfaceのblowupを全て含むdouble-point不等式

\(\epsilon:S\to S_0\) をminimal ruled modelとし、
底の種数をg、blowup回数をrとする。infinitely nearな中心も許す。
全変換exceptional divisorを \(E_i\) とすれば

\[
 H=\epsilon^*H_0-\sum_i m_iE_i,\qquad
 H_0^2=d+\sum_i m_i^2.
\]

\(k=H_0\cdot f\) を一般ruling fiberの次数とする。
Hのample性から \(k\ge1,m_i\ge1\)。
さらに

\[
 1\le m_i\le k-1. \tag{4}
\]

実際、各中心が属するruling fiberの全変換から
その時作られたexceptional divisorの全変換を引いたものは、非零有効divisorである。
Hとの交点が \(k-m_i>0\) なので (4) を得る。
これは中心がinfinitely nearでも、全変換のeffectivityにより同じである。

ruled surfaceのintersection formulaは

\[
 K_SH=-\frac{d+\sum m_i^2}{k}+2k(g-1)+\sum m_i,\quad
 K_S^2=8(1-g)-r,\quad \chi(\mathcal O_S)=1-g.
\]

これを(3)へ戻すと

\[
 0=d^2-10d+\frac{5d}{k}+(4-10k)(g-1)
       +\sum_i\left(\frac{5m_i^2}{k}-5m_i+2\right). \tag{5}
\]

この恒等式は有限Laurent算術で再展開した。

\(k\ge2\) とすると、最初の三項は \(d\le6\) により

\[
 d^2-10d+5d/k\le-3d/2<0.
\]

第四項は \(g\ge1\) により非正。
各summandについては

\[
 \frac{5m^2}{k}-5m+2-\left(\frac5k-3\right)
      =\frac5k(m-1)(m-k+1)\le0
\]

であり、(4)と \(k\ge2\) から各summandは高々 \(-1/2<0\) となる。
従って(5)は成立せず、\(k\ge2\) は全て不可能。
無限個のkを標本チェックした議論ではない。

従って \(k=1\)、(4)から \(r=0\)。
(5)は

\[
 d(d-5)=6(g-1)
\]

となり、\(d\le6,g\ge1\) で可能なのは
\((d,g)=(5,1),(6,2)\) だけである。
この有限整数リストも完全に列挙した。

## 6. degree6・genus2はvery amplenessとRRに反する

\(S=\mathbb P_B(F^\vee)\)、\(H=\mathcal O(1)\) と書く。
Fはrank2、\(\deg F=d\) である。
\((d,g)=(6,2)\) で \(H^1(B,F)\ne0\) と仮定すると、
Serre双対性から非零射 \(F\to K_B\) がある。
そのimageはline bundle \(M\subset K_B\)、\(\deg M\le2\) であり、Fのline quotientとなる。
対応するsection \(B\hookrightarrow S\) へのHの制限はMなので、
Hがvery ampleならMもvery ampleでなければならない。
genus2曲線の次数2以下のline bundleはvery ampleでなく、矛盾する。
よって \(H^1(B,F)=0\)、RRは

\[
 h^0(S,H)=h^0(B,F)=6+2(1-2)=4
\]

を与え、非退化 \(\mathbb P^4\) embeddingの五つの独立sectionと矛盾する。
残るのは \((d,g)=(5,1)\)、すなわち次数5の楕円scrollである。補題が証明された。□

## 7. very ampleな次数5scrollのbundleはstable

Bは楕円曲線、Fはrank2・次数5である。
Fがstableでなければ、飽和line subbundleの次数は3以上であり、
line quotient Mの次数は2以下となる。
これもHのsectionへの制限がvery ampleであることに反する。
従ってFはstableであり、
\(H^1(F)=\operatorname{Hom}(F,\mathcal O_B)^*=0\)、\(h^0(F)=5\) である。
元のembeddingは完全線形系による。

また

\[
 \operatorname{Hom}(\operatorname{Sym}^2F,\mathcal O_B)
   \hookrightarrow\operatorname{Hom}(F,F^*)=0.
\]

最後の消滅はstable bundleのslope \(5/2>-5/2\) から従う。
Serre双対性により \(H^1(\operatorname{Sym}^2F)=0\) となる。
従って

\[
 h^0(S,2H)=\deg(\operatorname{Sym}^2F)=15. \tag{6}
\]

finite étale splittingやSymmetric-power semistability定理を必要としない。

## 8. no quadric、3-regularityと全parameterのresolution

この滑らかな次数5楕円scrollはquadricに入らない。
rank5のquadric threefoldのPicard群はhyperplane classで生成されるため、
そこにある曲面の次数は偶数となる。
rank4の場合、一般ruling線は頂点を避け、一意の最大平面に入る。
それらの平面は二つの \(\mathbb P^1\) 族をなすので、
楕円ruling底から一方の \(\mathbb P^1\) へ写像を得る。
非定数なら次数2以上のため二ruling線が同一平面で交叉し、
定数ならscroll全体が一平面に入る。いずれもembeddingに反する。
rank3のquadricにも一つの \(\mathbb P^1\) 族の最大平面があり、同じ議論が適用される。
rank2以下は可約なので、既約・非退化なscrollを含めない。

よって \(H^0(I_S(2))=0\)、(6)と \(h^0(\mathbb P^4,\mathcal O(2))=15\) から
\(H^1(I_S(2))=0\) である。また

\[
 H^2(I_S(1))=H^1(F)=0,\qquad
 H^3(I_S)=H^2(\mathcal O_S)=0,\qquad H^4(I_S(-1))=0.
\]

従って \(I_S\) は3-regularであり、三次式で生成される。
regularityと一次完全性からprojective normalityを得る。
さらにRRとregularityから、全 \(n>0\) に対して
\(h^0(S,nH)=5n(n+1)/2\)。
座標環のHilbert seriesは

\[
 \frac{1-5t^3+5t^4-t^5}{(1-t)^5}.
\]

\(q(S)=1\) とprojective normalityから座標環のdepthは2であり、
最小ideal分解は長さ2。
3-regularityと最小性、Hilbert numeratorから

\[
 0\to\mathcal O(-5)\to\mathcal O(-4)^5
      \to\mathcal O(-3)^5\to I_S\to0. \tag{7}
\]

これは全stable rank2・次数5bundleに対する幾何学的証明であり、
generic有理函数体計算から特殊parameterへの飛躍ではない。

## 9. scalar q52認識

入口noteのBott公式と(7)から

\[
 h^0\bigl(\mathbb P^4,\Lambda^2T_{\mathbb P^4}\otimes I_S\bigr)
       \le h^2(\Omega^2_{\mathbb P^4})=1.
\]

Markarian--Polishchuk, arXiv:2310.18759v2, §3.2,
Lemmas 3.2.1,3.2.3の任意の楕円曲線とrank2・次数5 stable bundleに対する構成は、
同じtautological embedded scroll上で零になる非零 \(q_{5,2}\) を与える。
非零性は、その零locusが真の曲面であるという同lemmaからも分かる。
従って上の切断空間はちょうど一次元。
今回のpiもS全体で零なので、scalar \(q_{5,2}\) に一致する。
これで認識定理が証明された。□

## 10. この認識で捨てていない補集合

既約非正規五次の全支持をこの定理が尽くしたとはしない。
以下は明示的に残る。

* smooth Sでも一般rank2である非正規曲面。
* 特異な非正規曲面。小次数の非ACM曲面には入口noteのsyzygy必要条件だけを保持する。
* generic ordinary doubleでない横断型。
* smooth rank0・generic nodeでも、strict transformとEの交叉にvertical divisorがある型。

特にvertical成分が存在すると、その一般点は今回のhorizontal node計算の対象ではなく、
betaのpoleを消す根拠がない。この場合にq(S)正やscroll認識を主張しない。
同様に、任意の正規化や導手を滑らかと読み替えたわけではない。
今回除いた仮定と、なお残したclean条件を区別する。

有限算術は約0.0016秒、Python最大RSS13,312 KiBで完了した。
旧entrance・ordinary-doubleの有限関数とM2入力関数はAST比較で不変を確認し、
旧ordinary-double全算術の再実行結果も完全一致した。
旧immutable reportの埋込みsource bytesは保持し、互換mode追加理由を旧jobに明記した。
