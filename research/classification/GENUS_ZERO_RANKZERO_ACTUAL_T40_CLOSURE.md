# 次数五・種数零 rank 零全源の T40 帰属：特異 genus 一 base の部分正規化

2026-11-08、PRIMARY E。純代数幾何 producer proof。独立レビュー前。
新しい行列、GB、minor、jet、選択点又は normal form の計算は行わない。

**定理。** actual 被約既約 quintic \(D_0\subset\mathbb P^4\) を持つ rank 四 Poisson tensor
\(\pi_0\) の degree 五 rank 零支持が、正規化 sectional genus 零分岐に入るなら、
\([\pi_0]\) は閉成分 T40 に入る。
全 \(S(1,1,3)\)／\(S(1,2,2)\)、全五 section の射影、plane 導手の node／cusp、全 frame を保持する。
一般射影の閉包又は vector bundle の閉包だけによる主張ではなく、
全中心 tensor を一致させる実 ambient Poisson family を構成する。

## 1. 全 source と使用する旧定理

`RANK_ZERO_QUINTIC_SURFACE_GENUS_ZERO_NONDOMINANCE.md` とその
`RANK_ZERO_QUINTIC_SURFACE_GENUS_ZERO_REVIEW.md` の全源還元を使う。
normalization cone、degree 二 plane 導手 \(Q\subset S(1,1,3)\)、generic cusp \(S\) は既に空。
残る正規化は
\[
 Y_0=\mathbb P_{\mathbb P^1}(\widetilde E),\qquad
 \widetilde E=O(1)^2\oplus O(3)\quad\hbox{又は}\quad O(1)\oplus O(2)^2.
\]
quotient convention で \(L=O_{Y_0}(1)\)、\(L^3=5\)、\(h^0(L)=8\) である。
元の有限双有理射 \(f_0:Y_0\to D_0\) は任意の五次元部分空間
\(V\subset H^0(\widetilde E)\) とその指定された target frame による。
唯一の残余非正規 plane \(P\) の導手は
\[
 C_P=F_p+F_q\ (p\ne q),\quad\hbox{又は}\quad C_P=2F_p.
\]
以下、同じ記号 \(c=p+q\) 又は \(2p\) を \(\mathbb P^1\) 上の degree 二 divisor に使う。
全場合で
\[
 W=H^0(\widetilde E(-c))\subset V,\qquad \dim W=2,
 \qquad H^1(\widetilde E(-c))=0.                         \tag{1}
\]
plane の二つの linear equations がこの全 \(W\) を尽くす。
\(S\) は generic node で、正規化 tensor の全 source は
\[
 H^0(Y_0,\Lambda^2T_{Y_0}(-C_S))
 \simeq H^0(Y_0,\Omega^1_{Y_0}(C_P))=\mathbb C.          \tag{2}
\]
この一形式は base から来るので、全 tensor は \(Y_0/\mathbb P^1\) に vertical である。
plane 上の rank 零／二の両方、特別な射影中心、導手の codimension 二境界を追加仮定で除かない。

## 2. 五 section が nodal／cuspidal base 上の局所自由 bundle を定義する

### node の全 endpoint gluing

\(p,q\) を同一視した integral nodal plane cubic を \(B_0\)、その正規化を
\(n:\mathbb P^1\to B_0\) とする。両 fiber \(F_p,F_q\) は同じ plane \(P\) へ写る。
五 section の fiber evaluation は全射で、その kernel は両方とも \(W\)。
従って
\[
 V/W\xrightarrow{\sim}\widetilde E_p,
 \qquad V/W\xrightarrow{\sim}\widetilde E_q
\]
が fiber identification \(\varphi:\widetilde E_p\to\widetilde E_q\) を定める。
この \(\varphi\) で \(\widetilde E\) を貼り合わせた rank 三 bundle を \(E_0\) on \(B_0\) とする。
局所 node 環 \(A=\{(a,b)\in\mathbb C[[x]]\oplus\mathbb C[[y]]:a(0)=b(0)\}\) 上では
module は \(\{(u,v):v(0)=\varphi u(0)\}\simeq A^3\)。よって局所自由である。

(1) により二 fiber への全 evaluation \(H^0(\widetilde E)\to\widetilde E_p\oplus\widetilde E_q\) は全射。
貼合せ条件は独立三条件であり、\(h^0(E_0)=8-3=5\)、その全 section 空間はちょうど \(V\)。
また \(h^1(E_0)=0\)、\(n^*E_0=\widetilde E\)、\(\deg E_0=5\) である。
任意の \(\varphi\) を許し、cyclic 又は semisimple という条件は課さない。

### cusp の全 first-jet gluing

今度は \(2p\) を一点へ pinch した integral cuspidal plane cubic を \(B_0\) とする。
(1) により \(H^0(\widetilde E)\to H^0(\widetilde E|_{2p})\) は全射。
\(V/W\) の像は dimension 三で、fiber \(\widetilde E_p\) への評価は同型である。
局所 parameter \(z\) と bundle frame を選ぶと、この三次元像は
\[
 u(z)=u_0+zN u_0\pmod {z^2}
\]
という任意の \(3\times3\) matrix \(N\) の graph となる。
cusp 環 \(A=\mathbb C[[z^2,z^3]]\subset\mathbb C[[z]]\) 上の対応 module は
\[
 (I+zN)A^3\subset\mathbb C[[z]]^3,
\]
従って自由 rank 三である。
これは conductor square \(2p\to\{\mathrm{cusp}\}\) での代数的 bundle gluing であり、
formal frame は局所自由性を確認するためだけに使った。
有限 subscheme 上の上記 graph により代数的 bundle \(E_0\) on \(B_0\) が得られる。
同じ独立三条件の計算から \(h^0(E_0)=5\)、\(H^0(E_0)=V\)、\(h^1(E_0)=0\)、\(\deg E_0=5\)。
全 \(N\) を許し、Jordan 型、trace、rank を固定しない。

両場合とも、\(V\) は特異点でも \(E_0\) を生成する。
node では共通 fiber 評価、cusp では graph の constant term が同型だからである。
他の点では元の basepointfree 射に等しい。
また \(O_{\mathbb P(E_0)}(1)\) は finite surjective な正規化からの pullback が ample な \(L\) なので ample。
特に \(H^0(E_0^*)=0\)。後者は \(n^*E_0^*=\widetilde E^*\) が section を持たないことからも直接従う。

## 3. plane を先に貼った Gorenstein 部分正規化

\[
 Z_0=\mathbb P_{B_0}(E_0)
\]
と置くと、\(Y_0\to Z_0\) は有限正規化であり、五 section の射は
\(Y_0\to Z_0\xrightarrow{g_0}D_0\) と factor する。
\(g_0\) は ample な tautological bundle によるので有限、合成が双有理なので自身も双有理である。
\(B_0\) は integral Gorenstein genus 一 curve、\(\omega_{B_0}=O_{B_0}\)。
従って \(Z_0\) も integral Gorenstein で、smooth な \(\mathbb P^2\)-bundle over \(B_0\) である。
\(Z_0\) 自身が smooth 又は normal だとは主張しない。
これは plane の二 branch 又は cusp first jet を**内部で正しく貼った**部分正規化である。
\(\sigma/k\) が正規化上で正則という理由だけで plane を越える descent を仮定していない。

## 4. 全五 section の quotient を代数的に smooth する

\(B_0\) を plane cubic smoothing \(\mathcal B\to T\) の中心に置く。
例えば node／cusp の平面三次方程式に、特異点で非零の一般三次式の \(t\) 倍を加える。
零を残して縮めれば \(\mathcal B\) は regular、全 family は flat Gorenstein、
\(t\ne0\) の \(B_t\) は smooth elliptic curve、adjunction により \(\omega_{\mathcal B/T}=O_{\mathcal B}\)。

bundle だけでなく、指定された**全五 section の quotient**
\[
 0\longrightarrow Q_0\longrightarrow O_{B_0}^{,5}
    \longrightarrow E_0\longrightarrow0                       \tag{3}
\]
を変形する。\(Q_0\) は locally free rank 二である。
relative Quot の、この quotient を固定した障害は
\[
 \operatorname{Ext}^1(Q_0,E_0)=H^1(Q_0^*\otimes E_0)
  =\operatorname{Hom}(E_0,Q_0)^*=0.                            \tag{4}
\]
等号は Gorenstein curve \(B_0\) の Serre duality と \(\omega_{B_0}=O\) を使う。
両入力が locally free なので、特異点での非完全 module に duality を誤適用していない。
零性は \(Q_0\subset O^5\) と \(H^0(E_0^*)=0\) から従う。

(4) により relative Quot はこの点で smooth over \(T\)。
この標準 deformation statement は、局所 quotient の lifts を overlap 上で貼る差が
\(\mathcal Hom(Q_0,E_0)\) の Čech 一 cocycle となり、その \(H^1\) が零であることでも確認出来る。
Quot は finite presentation なので infinitesimal smoothness は actual smoothness を与える。
smooth morphism の局所 étale coordinates で相対座標を中心値に固定すれば、
零を含む étale-local **代数的曲線**上の quotient
\(O_{\mathcal B}^{,5}\twoheadrightarrow\mathcal E\) が得られる。
有限 jet 又は未代数化の formal 変形に止まっていない。

properness を使って基底を縮め、\(\mathcal E\) は全域 locally free、tautological bundle は相対 ample と出来る。
\(h^1(E_0)=0\) と base change より \(p_*\mathcal E\) は rank 五、指定の五 section はその frame となる。
各 \(E_t\) は degree 五 rank 三、ample／globally generated。
中心の任意の \(V\) とその target frame は (3) として固定したので失われない。

## 5. 全 fiber の実 quintic image と finite conductor

\(Z=\mathbb P_{\mathcal B}(\mathcal E)\) の五 section による射を
\(g:Z\to\mathbb P^4\times T\) とする。
全 fiber で tautological bundle は ample／globally generated なので有限射であり、\(g\) 自身も proper quasi-finite、従って有限。
各 image は非退化三次元で、\(L_t^3=5\) より
\(5=\deg(g_t)\deg(g_t(Z_t))\)。五の素数性と非退化性から全 fiber は双有理な integral quintic image。

total image \(D\) は regular ambient 中の horizontal integral prime divisor、従って Cartier。
その相対 degree は五。base line bundle を零の近くで trivialize すれば一本の degree 五 equation \(F_t\) を取れる。
horizontal prime は全 fiber を成分に持たないので \(F_t\) は零多項式にならない。
各 \(F_t\) は同 degree の integral quintic image を含むため、scheme としてそれに等しい。
従って全 actual image の flat hypersurface family を得て、\(\omega_{D/T}=O_D\)。

\(Z/T\) は Gorenstein／CM であり、total \(Z\) は regular な \(\mathcal B\) 上の projective bundle。
有限同次元射の \(B=g_*O_Z\) は \(A=O_D\) 上 MCM。
有限双対性は
\[
 \operatorname{Hom}_A(B,A)=g_*\omega_{Z/T},\qquad
 \operatorname{Ext}^{i}_A(B,A)=0\ (i>0).                         \tag{5}
\]
\(Z_0\) が非正規でもこれは有効であり、normality は MCM の代用として使っていない。
有限双有理な共通 fraction field 内で \(\varphi\mapsto\varphi(1)\) は左辺を
\(\mathfrak c_g=\{b\in B:bB\subset A\}\subset A\subset B\) と同定する。
従って自然な map \(\omega_{Z/T}\to O_Z\) を得る。
\(\operatorname{Ext}^1=0\) と \(0\to A\xrightarrow{t-t_0}A\to A_{t_0}\to0\) から
\[
 \operatorname{Hom}_A(B,A)/(t-t_0)
  =\operatorname{Hom}_{A_{t_0}}(B_{t_0},A_{t_0})                 \tag{6}
\]
で、evaluation at one も可換。中心を含め**正確な fiber conductor map**に特殊化する。
これは \(Y_0\) の full normalization conductor ではなく、\(Z_0\to D_0\) の partial conductor である。

## 6. plane descent を含む全 tensor の構成

\(Z\to\mathcal B\) は smooth relative dimension 二である。
\(\omega_{\mathcal B/T}=O\) の相対 dualizing generator を固定すると、
\(\omega_{Z/T}^{-1}\simeq\omega_{Z/\mathcal B}^{-1}=\Lambda^2T_{Z/\mathcal B}\)。
(5) の conductor inclusion をこの同型で読むと、全 family の**vertical bivector** \(P_Z\) が得られる。
base fiber が node／cusp でも \(T_{Z/\mathcal B}\) は rank 二 bundle である。
任意の vertical bivector は vertical 三次 exterior が零なので Poisson biderivation となる。
特異 base の dualizing differential を Kähler differential と同一視したのではない。

\(dg(P_Z)\) の係数は全て \(\mathfrak c_g\) に入るため、
\(\mathfrak c_gB\subset A\) により \(\Lambda^2T_{\mathbb P^4}|_D\) の正則 section \(\xi\) に descent する。
plane の node matching／cusp jet matching は \(Z_0\) の構成で既に満たされているので、
正規化上だけの正則性から gluing を外挿する隙間は無い。

中心の \((P_Z)_0\) を \(Y_0\) に pullback すると、relative tangent は base change で同型なので正則である。
\(S\) の generic preimage では \(Z_0=Y_0\) で、\(g_0\) の conductor は \(C_S\) に沿って零。
従って pullback は (2) の非零線に入る。
finite birational conductor map は generic に非零なので tensor は零ではない。
相対 dualizing generator の一つの定数を調整すれば、中心は指定された**元の全 source tensor** \(\pi_{Y_0}\) に等しい。

quintic restriction の exact sequence と Bott の消滅
\[
 H^0(\Omega^2_{\mathbb P^4})=H^1(\Omega^2_{\mathbb P^4})=0
\]
により、\(\xi\) は ambient bivector \(\Pi\) へ一意に正則 lift する。
中心は \(\pi_0\)。Jacobi は dense smooth \(D\) 上で零で、
\(H^0(\Lambda^3T(-5))=H^0(\Omega^1)=0\) より全 ambient で零である。
actual Pfaffian は \(c(t)F_t\)、\(c(0)\ne0\) なので、零を残して縮めると全 family が実 rank 四／被約既約 quintic となる。
全 section／conductor／ambient Poisson を lift した実曲線であり、projection 閉包だけの議論ではない。

## 7. smooth elliptic 全 fiber と T40 の閉性

ここで使う既受理証明の bundle-wide な部分を明示する。
smooth elliptic \(B\) 上の ample／GG rank 三 degree 五 bundle \(E\) が完全五 section により
quintic へ有限双有理に写り、その**全 conductor に沿って零**となる正規化 tensor が actual rank 四に lift するとする。
この条件なら、その tensor は T40 に入る。

stable \(E\) については、旧 genus 一認識の polarized-normalization 同型と tensor 一意性がそのまま bundle-wide に有効。
unstable \(E\) の全分類 \(E=A_2\oplus F_{2,3}\) は ample／GG のみを使う。
`UNSTABLE_ELLIPTIC_RANKZERO_ACTUAL_T40_CLOSURE.md` §§2–6 の逆 extension／全 conductor／ambient lift 証明も、
入口以降はこの bundle と \(H^0(\Omega^1)=\mathbb C\) の線だけを使う。
従ってこの lemma は conductor が一つの integral degree 五 surface かどうかには依存しない。
一般 fiber の conductor が split しないとの追加仮定を置いていない。

本構成の \(t\ne0\) は smooth elliptic \(B_t\)、ample／GG rank 三 degree 五 \(E_t\) であり、
conductor inclusion から作った tensor なのでこの全条件を満たす。
従って \([\Pi_t]\in\mathrm{T40}\)。T40 は閉であり、中心 \([\pi_0]\) も T40 に入る。
これが任意の node／cusp gluing、任意の \(V\)、任意の元の frame に適用するため、定理の全 source を覆う。

## 8. 113 型についての独立した弱い排除

実際、113 型は rank 四では空まで言える。
\(H^0(\widetilde E(-2))=H^0(O(1))\) は O(3) summand だけから来るので、
plane の二つの ambient linear equations \(\ell_0,\ell_1\) は pullback で
\(k(s,t)z\,s,k(s,t)z\,t\) となる。
(2) の tensor は vertical なので、\(\ell_0/\ell_1=s/t\) は \(D_0\) 上の非定数 Casimir。
projective field \(\ell_1\pi_0^\#d\ell_0-\ell_0\pi_0^\#d\ell_1\in H^0(T_{\mathbb P^4}(2))\) は \(D_0\) 上零となり、
\(H^0(T_{\mathbb P^4}(-3))=0\) から ambient でも零。
これは rank 四 symplectic 一般点に非定数 rational Casimir を作り矛盾する。
この排除は第2–7節の T40 包含に必要ではないが、全113補集合を単独で確認する。

## 9. 結論と保存範囲

rank 零 degree 五・normalization genus 零の全 rank 四 source は T40、113 は実際には空。
genus 一の既済 canonical／零 summand／cusp 排除と stable／unstable T40 を合わせると、
degree 五 rank 零導手 surface の全 rank 四分岐は T40 に入る。
新 component の増分は零。rank 二 restriction、他 degree／factor 型、全 E、総 component 数は未主張。
T36／T38 との未確認同定や、cyclic nodal normal form の Poisson descent を入力に使わない。
旧 source／report／主稿／global／Git は変更せず、同名 manifest に全依存を固定する。
