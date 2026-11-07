# 次数五 rank 零・正規化種数一の全認識：別作者による幾何監査

2026-11-07。producer は `generic_exhaustion`、reviewer は `log_manuscript`。
`RANK_ZERO_QUINTIC_SURFACE_GENUS_ONE_RECOGNITION.md`
（SHA256 `5dbc8052d4245b72d45904a0f92d479ce4e0c26a577ca2391a1f662bce35b768`）を
主義務 G の全源命題として独立に確認した。
以下の命題を受理する。CAS、数値標本、有限 jet 又は dual-span の target 同一視を使わない。
監査契約は `jobs/rank_zero_quintic_genus_one_review.json`、構造化記録は
`reports/CLASSIFICATION_RANK_ZERO_QUINTIC_GENUS_ONE_REVIEW_20261107.json` である。

**受理する命題。** \(D\subset\mathbb P^4\) が rank 四 Poisson tensor の実既約被約 quintic、
\(S\) が非正規支持の次数五 integral 成分、\(\pi|_S=0\) とする。
\(D\) の正規化 \((X,L)\) の sectional genus が一なら、全構造は T40 の閉像又は
全 PGL 飽和射影像の次元が24以下の族に含まれる。
\(S,X\) の smoothness、導手の normality、generic node、clean 又は vertical 不存在は仮定しない。

## 1. 新しい一次入力の実際の範囲

次の一次資料を reviewer が直接読解した。

* Y. Fukuma, *A lower bound for the second sectional geometric genus of quasi-polarized manifolds and its applications*,
  [2011 author/journal PDF](https://seminariomatematico.polito.it/rendiconti/69-1/73.pdf),
  Definition 1（印刷75頁）、Proposition 3(iii) と証明（76頁）。
  対象は **smooth quasi-polarized manifold** である。
  sectional genus 一なら、polarization の共通引戻しを持つ birational model は
  Gorenstein・Q-factorial・terminal の \(K+2H\sim0\) 型又は elliptic \(\mathbb P^2\)-scroll 型に尽くされる。
  singular \(X\) へ同定理を直接適用せず、smooth resolution に適用する。
* Greb--Kebekus--Kovács--Peternell,
  [arXiv:1003.2913v4](https://arxiv.org/pdf/1003.2913v4), Theorem 1.4（印刷3頁）。
  complex quasi-projective **klt** variety と log resolution に対し、regular \(p\)-form の押し出しは reflexive。
  よって \(X_{\rm reg}\) 上の一形式は resolution へ regular に延長する。
  arbitrary normal variety 又は単なる log canonical variety への正則延長とは読まない。
* M. F. Atiyah, *Vector bundles over an elliptic curve*,
  [1957 primary PDF](https://math.berkeley.edu/~nadler/atiyah.classification.pdf),
  Theorem 7 と Corollary（印刷434頁）。
  rank と degree が互いに素なら、indecomposable bundle は determinant により一意。
  ここで最後に使う rank は三、degree は五である。
  併せて Theorem 5（432頁）とその証明を読み、degree 零 indecomposable rank 二の自己 extension 型も照合した。

既知 \(q_{5,2}\) の幾何は Markarian--Polishchuk,
[arXiv:2310.18759v2](https://arxiv.org/html/2310.18759v2) §3.2 の冒頭、
Lemmas 3.2.1、3.2.3、3.2.4、Corollary 3.2.5、及び §3.4 の Lemma 3.4.1(i) と証明を直接読んだ。
以下 §6 では、`IRREDUCIBLE_NONNORMAL_QUINTIC_ENTRANCE.md` §4 の
**全 fiber に対する有限性の導出**をこの一次内容と再照合する。
一次資料に存在しない全 Poisson 分類を引用したのではない。

## 2. resolution からの全 ample model の二分岐

length 式から \(\delta_S=1\)、他に二次元非正規支持は無い。
\(\nu:X\to D\)、\(L=\nu^*O_D(1)\) は normal threefold と ample・globally generated Cartier divisor で、
\(L^3=5\)、\(K_X+C_X\sim0\) である。
二つの一般 hyperplane を引き戻した曲線は \(\operatorname{Sing}X\) を避ける滑らかな正規化曲線で、
その genus が一であることは既受理の Bertini・Hilbert 二階差分の議論による。
resolution \(Y\to X\) は \(X_{\rm reg}\) 上で同型に選べる。
同じ一般曲線は exceptional image を避けるので、\((Y,\rho^*L)\) の sectional genus も一である。

Fukuma をこの smooth quasi-polarized pair に適用する。
共通 smooth resolution \(Z\) 上で \(g^*L=h^*H\) が成り立ち、
\(X,V\) が normal なので \(g_*O_Z=O_X\)、\(h_*O_Z=O_V\)。
従って全 \(m\ge0\) について section ring は等しい。
\(H\) は globally generated でもある。実際、任意の \(v\in V\) の上の点を選び、
basepoint-free な共通引戻しの非零 section を取れば、対応する \(H\)-section は \(v\) で非零である。
よって section ring の ample model は元の \(X\) であり、morphism
\(f:V\to X\)、\(H=f^*L\) を得る。
これは birational model \(V\) を元の \(X\) と無条件で同一視する操作ではない。

## 3. canonical 型の全排除と closed inverse-residue

第一型では \(K_V\sim-2H\) を birational pushforward し \(K_X\sim-2L\) を得る。
従って \(K_X\) は Cartier である。
compatible canonical divisor を取ると \(K_V-f^*K_X\) は exceptional で numerically trivial。
negativity lemma を両符号に用いて零となる。
従って \(f\) は crepant、terminal \(V\) の resolution と比較して \(X\) は canonical、特に klt である。
\(X\) が Q-factorial 又は terminal そのものとは主張しない。

\(-K_V=2H\) は nef and big、\(V\) は terminal なので Kawamata--Viehweg により
\(H^1(O_V)=0\)。terminal の rational singularity 性と smooth common resolution の
birational invariance から \(H^1(O_Z)=0\)、従って \(H^0(\Omega_Z^1)=0\)。
ここで初めて、上で検査した klt 仮定の下で GKKP を使い
\[
 H^0(X,\Omega_X^{[1]})=H^0(Z,\Omega_Z^1)=0
\]
を得る。canonical Fano の未読個別分類又は twisted-reflexive 消滅に置き換えていない。

\(\delta=1\) の一般横断型は node 又は cusp に尽くされる。
node の全 branch quotient 計算は
`RANK_ZERO_SEXTIC_CONDUCTOR_SEPARATE_REVIEW.md` §5 と同じで、
rank 零から rational \(\pi_X\) は reduced conductor に沿って零となる。
cusp では \(A=K[[t^2,t^3]]\subset B=K[[t]]\)、\(C_X=2\Gamma\) である。
rank 零から base 同士の bracket は \((t^2,t^3)\) に入り、mixed 成分は
\(\{a,t\}=\{a,t^2\}/(2t)\in tB\)。
二つの base 方向と normal parameter の全成分を保持したまま、\(\pi_X\) は少なくとも一次で零となる。
係数体を Poisson 中心とは仮定しない。

四次元で \(\pi^{-1}\) を作る inverse-Pfaffian 恒等式は、reduced な actual \(D\) に沿って
高々一位極を持つ closed symplectic 二形式を与える。
smooth \(D\) 上の residue \(\alpha_D\) は closed であり、
dualizing generator \(\Omega_D\) との関係は
\(\alpha_D=c\,i_{\pi_D}\Omega_D\)、\(c\ne0\) である。
\(\omega_D=O_D\) と有限双対性から \(\Omega_D\) の rational pullback の divisor は正確に \(-C_X\)。
よって \(\alpha\) は node では divisorial pole 無し、cusp では唯一の \(\Gamma\) に高々 simple pole を持つ。

closed な simple-pole 一形式の \(dt/t\) residue は \(\Gamma\) 上で constant \(\lambda\) となる。
閉性は tangential \(da/t\) 項も排除する。
一般の二 \(|L|\) 切断曲線は \(\operatorname{Sing}X\) と codimension 二 bad locus を避け、
\(\Gamma L^2=5\) 個の交点を transverse に持つ。
その compact curve 上で residue theorem は \(5\lambda=0\) を与え、\(\lambda=0\) となる。
normal \(X\) の全 codimension 一で pole が消えるので \(\alpha\in H^0(\Omega_X^{[1]})\)。
しかし actual reduced \(D\) の一般点で \(\pi_D\) は非零 rank 二であり、\(\alpha\ne0\)。
上の GKKP 消滅と矛盾する。singular \(S\) の reflexive 形式を resolution へ無条件に延長したのではない。

## 4. elliptic nef/GG source の零 summand と全 unstable 型

第二型は \(V=\mathbb P_B(E)\)、\(H=O(1)\) で、\(B\) は smooth elliptic curve、
\(E\) は globally generated rank 三 degree 五である。
非零 \(E\to O_B\) があれば、その image \(O_B(-Z)\) は globally generated なので \(Z=0\)。
global generation の評価写像を合成すると constant coefficient の \(O_B^N\to O_B\) が全射であり、
一つの global section が1へ写るためこの quotient は split する。
繰り返し \(E=O_B^r\oplus F\)、\(H^0(F^*)=0\) を得る。

この \(F\) は ample である。
globally generated な tautological morphism に contracted curve があれば、
その正規化 \(T\to B\) は surjective で、degree 零 globally generated line quotient は \(O_T\)。
global evaluation の \(H^0(F)\to H^0(O_T)\) は constant functional で、
全 \(T\) 上の evaluation kernel を消す。
\(T\to B\) の surjectivity から全 \(B\) の evaluation kernel も消し、\(F\to O_B\) が出来て矛盾する。
従って tautological morphism は finite、\(O(1)\) は ample となる。

\(r>0\) の場合、**完全 degree 一系**の image は \(O_B^r\) 方向を vertex とする cone である。
元の五 section はその image からの linear projection である。
\(\nu\) が finite かつ basepoint free なので、射影中心は vertex を含む complete image を避け、
vertex がそのまま \(D\) の vertex へ写る。
唯一の degree 五非正規成分 \(S\) も cone となり既受理の全源排除に反する。
\(r>0\) で complete section 数が五を超え得る点も、この射影の記述で保持した。
よって \(r=0\)、\(H\) は ample、\(f\) は finite birational なので \(X=V\) である。

smooth \(X=\mathbb P_B(E)\) では \(C_X=-K_X=3L-p^*\det E\) の \(L\) 係数は奇数。
cusp の \(C_X=2\Gamma\) は不可能で、全一般横断型は node となる。
従って全源は
\[
 H^0(X,\Lambda^2T_X(-C_X))=H^0(X,\Omega_X^1)=\mathbb C.
\]
ample \(E\) の dual に section は無いので \(h^1(E)=0\)、elliptic RR から \(h^0(E)=5\)。
元の \(\nu\) は complete 五 section の射である。
ambient tensor の restriction も
\(H^0(\mathbb P^4,\Lambda^2T(-5))=H^0(\mathbb P^4,\Omega^2)=0\)
により injective で、各 polarized bundle の projective tensor choice は高々一点。

ample/GG rank \(r\) quotient \(Q\) は \(h^1(Q)=0\)、\(h^0(Q)=\deg Q\ge r+1\) を満たす。
等号 \(h^0(Q)=r\) なら evaluation \(O_B^r\to Q\) が同型となり ample に反するからである。
従って rank 二 destabilizer の degree 四以上は quotient line degree 一以下を作るため存在しない。
rank 一 destabilizer は degree 二に限られ、rank 二 quotient は degree 三となる。
その quotient が unstable でも degree 一以下の line quotient を作るため、stable \(F_{2,3}\) である。
\[
 \operatorname{Ext}^1(F_{2,3},A_2)
 =H^0(F_{2,3}\otimes A_2^{-1})^*=0
\]
は stable degree \(-1\) bundle の section 零性から従い、全 unstable 型は
\(E=A_2\oplus F_{2,3}\) だけである。
逆向き Hom は零、\(h^0(F_{2,3}^*\otimes A_2)=1\) なので \(h^0(\operatorname{End}E)=3\)。
actual \(\operatorname{Aut}(E)/\mathbb C^*\) は二次元である。

\(B,A_2,F_{2,3}\) の parameter は各々一以下で全三以下。
stable \(F\) の局所 moduli 上界は \(h^1(\operatorname{End}F)=1\) からも確認出来る。
base translation は degree 五の determinant を一次元で動かし、その stabilizer は有限。
target frame 24、bundle automorphism 二と translation 一を実 fiber として数えると全像は
\(24+3-2-1=24\) 以下。
この作用の kernel は frame 付き complete finite birational map を固定する automorphism であり、
dense open で identity、従って全体で identity である。
tensor の stabilizer 零性を仮定した商次元ではない。

## 5. 既知 \(q_{5,2}\) の finite normalization を全点で再確認

任意の elliptic curve \(B_0\) と stable rank 二 degree 五 bundle \(W\) に対し、
MP §3.2 は \(\mathbb P H^0(W)^*\) 上の既知 tensor を扱う。
その zero support \(S_0=\mathbb P_{B_0}(W^\vee)\) は Lemma 3.2.3 により smooth embedded scroll。
degeneracy hypersurface \(D_0\) は \(A=\operatorname{Pic}^2(B_0)\) に沿う plane の和である。
構成された smooth \(\mathbb P^2\)-bundle \(f:X_0\to D_0\) は Lemma 3.2.4 により \(S_0\) 外で同型。
ここで止めず、\(S_0\) 上の全 fiber を確認する。

Lemma 3.2.1(i) の middle bundle は
\(O(p)\oplus W'\)、\(W'\) は rank 二 degree 四 semistable、\(\dim\operatorname{End}W'=2\) である。
fiber の点は degree 二 line bundle への splitting choice に対応する。
degree 二から \(O(p)\) への Hom は零なので、全 choice は \(W'\) に入る。
\(W'=L_1\oplus L_2\)、\(L_1\not\simeq L_2\) なら degree 二 subline はこの二つだけ。
indecomposable なら Atiyah の degree 零 rank 二を twist した非自明自己 extension で、
degree 二 subline は一つだけである。この自己 extension も MP Lemma 3.4.1 の証明に明記される。
従って全 fiber は有限集合であり、proper \(f\) は quasi-finite、ゆえに finite。
smooth \(X_0\) からの birational finite map なので **全 \(D_0\) の正規化**である。
特殊四点又は ramification を一般二点と取り違えていない。

\(L_0=f^*O_{D_0}(1)\) は finite pullback だから ample かつ GG。
quotient convention で \(X_0=\mathbb P_A(E_0)\)、\(L_0=O(1)\) とすれば rank は三、
\(\deg E_0=L_0^3=5\)、\(h^0(E_0)=5\) である。
元の target が complete 五 section の target であることがここで従う。
MP の零集合と off-scroll 同型から \(\operatorname{Sing}(D_0)_{\rm red}=S_0\)。
後半は、零点では Pfaffian の一次微分も零となることから従う。
MP の plane family は irreducible であり、\(S_0\) が \(\mathbb P^4\) を span するので
actual quintic は hyperplane の五重ではなく、既約被約 quintic である。

## 6. 安定 rank 三 bridge と全 T40 認識

上の \(E_0\) は stable である。
もし unstable なら §4 の完全分類から二次元の actual polarized automorphism があり、
complete finite birational map を通じて \(\operatorname{Aut}_{\rm PGL_5}(D_0)\) に faithful に写る。
しかし polarized \(S_0\) の automorphism group は零次元である。
Albanese ruling は保存され、固定 base 上の kernel は stable \(W\) の
\(\operatorname{Aut}(W)/\mathbb C^*=0\)。base translation は \(\det W\) を保存しなければならず、
degree 五だからその部分群は finite、origin-preserving base automorphism も finite である。
\(\operatorname{Aut}_{\rm PGL_5}(D_0)\) は singular support \(S_0\) を保存し、
\(S_0\) 上で identity なら irreducibility と nondegeneracy から ambient identity となる。
従ってこの群も零次元であり矛盾する。
rank 三 degree 五は互いに素なので strictly semistable の第三の場合も無い。

任意の elliptic \(B\) を \(\operatorname{Pic}^2(B_0)\) と同定出来るように \(B_0\) を選ぶ。
任意の stable rank 三 degree 五 \(E\) に対し、degree 五 determinant 上の translation は
\([5]:B\to B\) の isogeny なので、\(\det(t_b^*E_0)=\det E\) となる \(b\) が存在する。
Atiyah Theorem 7 の coprime Corollary を **二つの stable rank 三 bundle** に適用して
\(E\simeq t_b^*E_0\) を得る。
従って polarized normalization は同型であり、両方の complete 五 section の image は実 PGL 同値。
§4 の一意な projective \(\pi_X\) と ambient restriction injectivity から、tensor 自身も既知 T40 と一致する。

rank 二の \(W\) と rank 三の \(E\) を同じ bundle としたり、
dual-span だけから二つの \(\mathbb P^4\) target を同一視したりしていない。
全 fiber の finite normalization、polarized stabilizer、rank 三 determinant の一意性、
最後の complete-section map の順に実際の bridge を作った。

## 7. 最終判定と保存範囲

canonical 型は空、elliptic の零 summand 型は cone 排除、
全 unstable 型は像24以下、全 stable 型は T40 となる。
従って producer 命題を **degree 五 rank 零・正規化 sectional genus 一の全 G 分岐として受理**する。
genus 零は別の既受理全射影 bound、degree 六は別の全排除と合わせられるが、
plane、quadric cone、rank 二 restriction、全 G/E の完了又は成分数の変更はこの証明から主張しない。
保存済み146頁 TeX/PDF、global 台帳、Git は今回変更しない。
