# 滑らかな正規化と既約滑らかな導手を持つordinary-double五次の条件付き認識

2026-11-07。Eの既約非正規五次に対する条件付き命題である。
Fujita型の未確認の分類定理を使わず、導手上の二形式、有限群商、
楕円scrollの初等的なcohomologyからT40への認識を導く。
一般の既約非正規五次が以下の仮定を満たすとは主張しない。

有限算術契約は jobs/irreducible_nonnormal_ordinary_double.json、
実行器は code/classification_nonnormal_quintic_entrance.py の --ordinary-double、
証明書は reports/CLASSIFICATION_IRREDUCIBLE_NONNORMAL_ORDINARY_DOUBLE_RECOGNITION_20261107.json。
以下の幾何学的議論は計算証明書やgateのPASSとは別である。

## 1. 仮定と結論

複素射影四空間のrank4 Poisson構造を \(\pi\)、その実退化因子を
既約・被約・非正規五次 \(D\) とする。
\(\nu:X\to D\) は正規化、\(L_X=\nu^*\mathcal O_D(1)\) とする。
以下を仮定する。

1. \(X\) は滑らか。
2. 正規化上の導手divisor \(C\) は被約・滑らか・既約。
3. 非正規曲面の一般点で \(D\) は横断ordinary doubleである。
4. \(\pi^{-1}\) の \(D_{\rm reg}\) 上のresidueを \(X\) に移した有理閉一形式
   \(\alpha\) は、divisorialには \(C\) に沿う高々一位の極しか持たない。

仮定4は一般ordinary-double点の局所logarithmic表示から確認できるが、
高位極の可能性がある別の横断型へ無断で拡張しないため明記した。

**条件付き認識定理（今回証明、独立読解監査対象）。**
このとき \(C\) はabelian surfaceであり、導手image
\(S\subset D\subset\mathbb P^4\) は滑らかな次数5の楕円scrollである。
ある楕円曲線 \(B\) とrank2・次数5のstable bundle \(F\) により

\[
 (S,\mathcal O_S(1))\simeq
 (\mathbb P_B(F^\vee),\mathcal O_{\mathbb P(F^\vee)}(1)).
\]

その完全idealは

\[
 0\longrightarrow\mathcal O(-5)
 \longrightarrow\mathcal O(-4)^5
 \longrightarrow\mathcal O(-3)^5
 \longrightarrow I_S\longrightarrow0. \tag{1}
\]

従って \(\pi\) は、この埋込み楕円scrollに対応する
Feigin--Odesskii \(q_{5,2}\) のscalar倍であり、T40に属す。
genericな楕円parameterのM2計算を特殊parameterへ延長することは、この証明の入力ではない。

## 2. 二枝のresidueが導手のcotangent coframeを与える

一般ordinary-double点の局所方程式を \(uv=0\)、接方向を \(z,w\) とする。
\(\omega=\pi^{-1}\) は閉じており、SNC divisorに沿ってlogarithmicである。
log基底 \(du/u,dv/v,dz,dw\) で、そのskew行列を

\[
 W=\begin{pmatrix}
 0&\lambda&a&b\\
 -\lambda&0&c&d\\
 -a&-c&0&e\\
 -b&-d&-e&0
 \end{pmatrix},\qquad
 \operatorname{Pf}W=\lambda e-ad+bc \tag{2}
\]

と書く。\(\operatorname{Pf}W\) は一般二重点でunitである。
これは単に \(D\) 外で非退化というだけでなく、
\(D\) が実際の被約退化因子であることから従う。
実bivectorは

\[
 \pi=\operatorname{diag}(u,v,1,1)\,W^{-1}\,
             \operatorname{diag}(u,v,1,1).
\]

その \(u=v=0\) 上の接方向成分は
\(-\lambda(\operatorname{Pf}W)^{-1}\partial_z\wedge\partial_w\) である。
従ってこの二重曲面上のrank0は \(\lambda=0\) と同値である。
全inverse恒等式は有限算術でも分母を払って32成分検査した。

導手の二枝を交換する有理involutionを \(\iota\) と書く。
閉一形式 \(\alpha\) の \(C\) 上のresidueは定数であり、枝交換はその符号を反転させる。
\(C\) は既約なので、この定数は零である。
閉性と一位の極から \(\alpha\) はlogarithmicであり、
residue零なら極は消える。\(X\) が滑らかなのでcodimension二に隠れた有理形式の極もなく、
\(\alpha\in H^0(X,\Omega_X^1)\) となる。

\[
 \eta=\alpha|_C
\]

と置く。局所的に二枝上で得る接方向residueは
\(a\,dz+b\,dw\)、\(c\,dz+d\,dw\) である。
\(\lambda=0\) と (2) のunit性から

\[
 \eta\wedge\iota^*\eta\ne0. \tag{3}
\]

ここでinvolutionの正則延長と導手商の正規性は次節で証明する。
特に、この枝ではrank0は追加の無根拠な仮定ではなく、residueの符号から従う。
同時に、rank0をPoisson構造の不存在と取り違えていない。

有限双対性から \(\omega_X\simeq\mathcal O_X(-C)\) である。
随伴公式は \(K_C\sim0\) を与えるので、(3)は至る所非零のcanonical sectionとなる。
従って \(\Omega_C^1\) は二つの大域形式で自明化される。
よって \(c_2(C)=0\)、Noether公式から
\(\chi(\mathcal O_C)=0\) となる。
\(K_C\sim0\) より \(p_g(C)=1\)、従って \(q(C)=2\) である。
Albanese mapは至る所étaleで、同次元のabelian surfaceへの有限étale射となる。
原点を選べば \(C\) 自身がabelian surfaceである。

## 3. 導手imageは商の正規化だけでなく商そのもの

**導手gluing補題。** 上の仮定の下で

\[
 S=C/\iota
\]

がschemeとして成立する。特に \(S\) は正規である。

**証明。** 導手の定義から \(\mathcal O_S\hookrightarrow\nu_{C*}\mathcal O_C\) であり、
\(C\) が被約なので \(S\) も被約である。
\(C\to S\) は有限で一般次数2、\(C\) は正規である。
従って \(C\) は \(S\) の函数体の二次拡大内での整閉包であり、
そのGalois involutionは \(C\) 全体の正則involutionへ延びる。
\(T=C/\iota\) は \(S\) の正規化で、有限双有理射 \(T\to S\) を持つ。

導手で割った正規化完全列は

\[
 0\to\mathcal O_S\to\nu_{C*}\mathcal O_C
          \to M\to0,\qquad
 M=\nu_*\mathcal O_X/\mathcal O_D. \tag{4}
\]

\(D\) はCM三次元、正規化は \(S_2\) なので深さ補題により
\(M\) は純二次元である（入口noteの命題1）。
一方 \(\mathcal O_T/\mathcal O_S\) は (4) 右端 \(M\) の部分層であり、
有限双有理性からその支持は高々一次元である。
純性に反するのでこの部分層は零、従って \(T=S\) である。□

この段階を省き、商が正規化であることだけで導手imageを滑らかと呼ぶ必要はない。
一般ordinary-double型、normal conductor、\(D\) のCM性が全て使われている。

\(\eta,\iota^*\eta\) の独立性から
\(H^0(C,\Omega_C^1)\) 上の \(\iota\) の固有空間次元は \(1+1\) である。
固定点があればその接空間固有値も \(+1,-1\) であり、
固定locusは滑らかな曲線、局所商はreflection商で滑らかとなる。
固定点がなければ商射はétaleである。
従って両場合とも \(S\) は滑らかで

\[
 q(S)=\dim H^1(C,\mathcal O_C)^\iota=1. \tag{5}
\]

## 4. 曲面次数は5以下、埋込みは非退化

一般の二hyperplane切断の正規化種数を \(g\)、\(d=\deg S\) とする。
導手は既約、一般横断型はnodeなので入口noteの式は

\[
 d=6-g. \tag{6}
\]

また \(0\to K_X\to\mathcal O_X\to\mathcal O_C\to0\) のcohomologyとSerre双対性から

\[
 2=h^1(\mathcal O_C)\le h^1(\mathcal O_X)+h^2(K_X)=2q(X).
\]

従って \(q(X)\ge1\) である。
一般のample smooth切断 \(H\)、さらに \(\Gamma\) に対してKodaira消滅を使うと

\[
 H^1(\mathcal O_X)\simeq H^1(\mathcal O_H)
       \hookrightarrow H^1(\mathcal O_\Gamma).
\]

よって \(g\ge q(X)\ge1\)、従って

\[
 d\le5. \tag{7}
\]

\(S\) がhyperplane \(\mathbb P^3\) に入れば、面としてはそのhypersurfaceであり、
標準完全列から \(h^1(\mathcal O_S)=0\) となる。
さらに低いspanの場合も \(S=\mathbb P^2\) で同じである。(5)に反する。
従って \(S\subset\mathbb P^4\) は非退化であり、
\(\ell=\mathcal O_S(1)\) について

\[
 h^0(S,\ell)\ge5. \tag{8}
\]

## 5. involutionがfreeなら次数だけで矛盾

\(\iota\) がfreeなら \(C\to S\) はétale二重被覆である。
canonical characterは \(\det(\iota^*|\Omega_C^1)=-1\) なので、
\(K_S\) は非自明な2-torsion、\(p_g(S)=0\)、\(\chi(\mathcal O_S)=0\) となる。
ampleな \(\ell\) に対して \(\ell\otimes K_S^{-1}\) もampleなのでKodaira消滅とRRは

\[
 h^0(S,\ell)=\chi(S,\ell)=\frac{\ell^2}{2}=\frac d2.
\]

(7),(8)と両立しない。従って固定曲線が存在する。
bielliptic surfaceの分類を引用することなく、このfree商は排除された。

## 6. 固定曲線商は有限étale base changeで積になる

固定点を原点として \(C\) をabelian varietyと見れば、\(\iota\) はgroup automorphismである。
その \(+1,-1\) elliptic subvarietyを \(E_+,E_-\) とすると、additionはisogeny

\[
 E_+\times E_-\longrightarrow C
\]

となる。kernel \(K\) の各因子への射影は単射であり、
\(K\) は各因子の2-torsionに入る。
よって

\[
 S=(E_+\times(E_-/\{\pm1\}))/K
      =(E_+\times\mathbb P^1)/K.
\]

\(K\) は \(E_+\) 上で非零translationによりfreeに作用する。
従ってある楕円曲線 \(B=E_+/\operatorname{pr}_+K\) 上の \(\mathbb P^1\)-bundle
\(p:S\to B\) が得られ、
\(E_+\to B\) による有限étale base changeで \(E_+\times\mathbb P^1\) になる。

そのcoverの次数を \(N\) とする。ample \(\ell\) の引戻しは
\(M\boxtimes\mathcal O_{\mathbb P^1}(k)\) と書け、
\(m=\deg M>0,\ k\ge1\) である。従って

\[
 Nd=2mk,\qquad N(\ell K_S)=-2m.
\]

cover上で高次cohomologyが零であり、traceにより \(H^i(S,\ell)=0\) (\(i>0\)) である。
RRから

\[
 h^0(S,\ell)=\frac{d(k+1)}{2k}. \tag{9}
\]

(7),(8),(9)は \(d=5,k=1\) を強制する。
実際 \(d(k+1)\ge10k\) より \(k\le d/(10-d)\le1\)、その後 \(d\ge5\) である。

## 7. stable rank2次数5bundleと完全埋込み

\(F=p_*\ell\) と置くとrank2bundleで、相対評価により
\(S=\mathbb P_B(F^\vee)\)、\(\ell=\mathcal O(1)\) である。
上のétale coverで \(F\) の引戻しは \(M\oplus M\) になるので \(F\) はsemistable。
\(\deg F=\chi(B,F)=\chi(S,\ell)=5\) だから、rankとの互いに素性によりstableとなる。
さらに \(h^0(S,\ell)=5\) なので、元の \(\mathbb P^4\) embeddingは完全線形系による。

このstable bundleのtautological embedding自体も直接検査できる。
\(B\) 上の任意の長さ2のsubscheme \(Z\) に対して \(F(-Z)\) はstable、次数1であり、
Serre双対性から \(H^1(F(-Z))=0\) となる。
従って二fiberの値と一fiberでの一jetを分離し、fiber上の \(\mathcal O(1)\) と合わせて
tautological line bundleはvery ampleである。

## 8. このscrollはquadricに入らない

一般に、滑らかに埋め込まれた次数5楕円scrollはquadricに入らない。
これをprojective-normalityを先取りせず証明する。

rank5のquadric threefoldは滑らかでPicard群がhyperplane classで生成されるので、
その上の曲面の次数は偶数となり、次数5と矛盾する。

rank4の場合、頂点は一点、最大平面は二つの \(\mathbb P^1\) 族をなす。
scrollのruling線は互いに交わらないので、一般ruling線は頂点を避ける。
その線は一意の最大平面に入り、楕円ruling底から一方の \(\mathbb P^1\) へ写像を得る。
非定数なら次数は少なくとも2で、二つの異なるruling線が同じ平面に入り交わってしまう。
定数ならscroll全体が一平面に入る。いずれも矛盾する。

rank3の場合も、頂点線を含む最大平面が一つの \(\mathbb P^1\) 族をなす。
頂点線そのものでない各ruling線は一意の最大平面に入り、同じ議論が適用される。
rank2以下は可約quadricなので、既約で非退化なscrollを含めない。
従って \(H^0(I_S(2))=0\) である。

## 9. 全stable elliptic scrollのresolutionとscalar認識

上のétale coverを使うと、全 \(n>0\) に対して

\[
 H^{>0}(S,\ell^n)=0,\qquad
 h^0(S,\ell^n)=\frac{5n(n+1)}2.
\]

特に \(h^0(\ell^2)=15=h^0(\mathbb P^4,\mathcal O(2))\) であり、
前節のquadric不存在から \(H^1(I_S(2))=0\) となる。
また

\[
 H^2(I_S(1))=H^1(S,\ell)=0,\quad
 H^3(I_S)=H^2(S,\mathcal O_S)=0,\quad
 H^4(I_S(-1))=0.
\]

従って \(I_S\) は3-regularである。
一次・二次generatorはなく、三次generatorは \(35-30=5\) 個。
regularityから全次数でprojective normalityを得る。
座標環のHilbert seriesは

\[
 1+\frac{5t}{(1-t)^3}
   =\frac{1-5t^3+5t^4-t^5}{(1-t)^5}. \tag{10}
\]

さらに \(h^1(\mathcal O_S)=1\) とprojective normalityから座標環のdepthは2である。
従って \(I_S\) の最小自由分解の長さは2。
3-regularityと最小性から各項のshiftは3,4,5に限られ、
(10)はrankを5,5,1と決める。これで (1) が全楕円parameterに対して証明された。

\(E=\Lambda^2T_{\mathbb P^4}=\Omega^2(5)\) とする。
入口noteのBott/syzygy議論、または (1) の短完全列を直接使うと

\[
 h^0(EI_S)\le h^2(E(-5))=1. \tag{11}
\]

Markarian--Polishchuk, arXiv:2310.18759v2, Lemmas 3.2.1,3.2.3は、
任意の楕円曲線とrank2・次数5 stable bundleに対応する
\(q_{5,2}\) が、この同じtautological scroll上で零になることを与える。
ここでbundleの任意性は同論文§3.2冒頭のstable rank2・次数5という設定の範囲である。
Lemma 3.2.1(i)の零locusが真の曲面であることから、そのbivectorは非零である。
よって (11) の次元はちょうど1である。
今回の \(\pi\) は局所residue計算から \(S\) の一般点で零であり、
\(S\) が既約・被約で \(\pi\) が正則なので \(S\) 全体で零となる。
従って \(\pi\in H^0(EI_S)\) であり、scalar \(q_{5,2}\) に一致する。□

## 10. 残す補集合と検証状態

この定理は、滑らかな正規化、滑らかで既約な導手、一般ordinary-double、
正規化residueの一位極という仮定を持つ枝だけを認識する。
特異な正規化、特異または可約な導手、cuspidalその他の高位極型は残す。
それらが未知G優越成分の一般元になれないという次元評価は証明していない。

有限算術は全log-skew inverse、rank0/coframe恒等式、
RR不等式の全正整数解、Hilbert numeratorのみを確認する。
global gluing、abelian化、有限étale quotient、quadric排除、
regularityからのresolutionは本文の幾何学的証明であり、
計算PASSから自動的に従ったとは書かない。
