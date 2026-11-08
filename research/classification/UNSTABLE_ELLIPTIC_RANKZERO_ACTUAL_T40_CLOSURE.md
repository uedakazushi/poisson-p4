# 不安定 elliptic 正規化の全点は実 Poisson 曲線で T40 に入る

2026-11-08、PRIMARY E。純代数幾何による全 source の点ごとの包含。
新しい行列、GB、minor、jet 又は選択点の計算は無い。
旧 genus 一の G 証明の dimension 上界を変更せず、残っていた E の矢印を証明する。

**定理。** 一般 rank 四の Poisson tensor \(\pi_0\) on \(\mathbb P^4\) の実退化因子が
被約既約 quintic \(D_0\) であり、degree 五 rank 零導手曲面の genus 一分岐が
不安定 elliptic 正規化を持つとする。そのとき \([\pi_0]\) は既知の閉成分 T40 に入る。
各不安定 bundle と完全五 section の各 frame に対し、\(\pi_0\) を中心とする
**実際の ambient Poisson tensor の代数的曲線**を構成する。
bundle の安定型閉包だけからの推論ではない。

## 1. 受理済みの全入口と記法

`RANK_ZERO_QUINTIC_SURFACE_GENUS_ONE_RECOGNITION.md` §§3–6 と
`RANK_ZERO_QUINTIC_SURFACE_GENUS_ONE_SEPARATE_REVIEW.md` §§3–6 を使う。
canonical 分岐、零 summand と nef 非 ample 境界、cusp はこの入口で既に排除されている。
ここで \(B\) は smooth elliptic curve、
\[
 E_0=A\oplus F,\qquad \operatorname{rk}A=1,\ \deg A=2,
 \qquad \operatorname{rk}F=2,\ \deg F=3,
\]
であり、\(F\) は stable。\(A_2\) の添字二は degree で、rank ではない。
全不安定型はこれだけであり、逆向きの元の extension は split することも受理済みである。
\(E_0\) は ample／globally generated、\(h^0(E_0)=5\)、\(h^1(E_0)=0\)。

正規化 \(f_0:Y_0=\mathbb P_B(E_0)\to D_0\subset\mathbb P^4\) は
\(L_0=O_{Y_0}(1)\) の**完全**五 section による finite birational map である。
商 line の convention を使い、\(p_*O_Y(1)=E\)、\(L_0^3=\deg E_0=5\) とする。
divisorial conductor \(C_0\) は \(-K_{Y_0}\)、一般 transverse type は node であり、
rank 零から正規化された tensor は conductor に沿って零となる。
従って全 source は
\[
 \pi_{Y_0}\in H^0(Y_0,\Lambda^2T_{Y_0}(-C_0))
  \simeq H^0(Y_0,\Omega^1_{Y_0})=p^*H^0(B,\Omega^1_B)=\mathbb C.       \tag{1}
\]
これは一般点だけでなく、全 divisor に沿う零から smooth \(Y_0\) 上で得た global section である。
codimension 二の導手境界を捨てたとの仮定は無い。
ambient restriction は \(H^0(\Lambda^2T_{\mathbb P^4}(-5))=H^0(\Omega^2_{\mathbb P^4})=0\)
により injective なので、この一つの projective source を実際に lift すれば全点を覆う。

## 2. 逆向き extension は全非零 parameter で安定

今度は元と**逆向き**の extension を取る。
stable rank 二 degree \(-1\) bundle \(F\otimes A^{-1}\) は section を持たず、elliptic RR から
\[
 \operatorname{Ext}^1(A,F)=H^1(F\otimes A^{-1})\simeq\mathbb C.
\]
非零 class \(\eta\) を固定し、\(T=\mathbb A^1\) 上で
\[
 0\longrightarrow F\longrightarrow E_t\longrightarrow A\longrightarrow0,
 \qquad [E_t]=t\eta                                             \tag{2}
\]
を作る。transition matrix の off-diagonal cocycle を \(t\eta\) とすれば、
これは \(B\times T\) 上の algebraic vector bundle \(\mathcal E\) であり、特別 fiber は指定の \(E_0\)。

\(t\ne0\) なら \(E_t\) は stable である。実際、degree 二以上の line subbundle \(M\) があれば、
\(M\to A\) が零なら \(M\subset F\) となって \(F\) の stability に反する。
非零なら \(\deg M\le2\)、従って \(M\simeq A\) で quotient への写像は同型となり、(2) を split して矛盾する。
また degree 四以上の saturated rank 二 subbundle があれば、その quotient は degree 一以下の line \(Q\)。
\(\operatorname{Hom}(F,Q)=0\) なので \(E_t\to Q\) は \(A\to Q\) を通るが、これも degree に反する。
これが全 destabilizing rank 一／二の補集合を尽くす。

この証明は \(A,F\) の特殊な determinant 又は section を一切反転しない。
非零 \(t\) は全て stable であり、\(t=0\) が任意の元の不安定 source である。

## 3. 完全五 section と全 normalization map を lift する

以下、\(T\) を零を含む Zariski open に縮める。
\(h^1(E_0)=0\) と cohomology/base change より \(p_{T*}\mathcal E\) は rank 五 vector bundle で、
全 fiber の \(H^0(E_t)\) と一致する。
\(L_0\) の ampleness と global generation は開条件なので、
\(Y=\mathbb P_{B\times T}(\mathcal E)\) の \(L=O_Y(1)\) は全 fiber で ample／globally generated と出来る。
完全 section bundle をさらに零の近くで trivialize し、中心の任意の指定 frame をそのまま延長する。
従って全五 section による射
\[
 f:Y\longrightarrow\mathbb P^4\times T
\]
が得られる。有限個の section だけを選び直したり、元の frame を失ったりしない。

各 fiber の射は ample な \(L_t\) で曲線を縮約しないから有限である。
properness と fiberwise quasi-finiteness により \(f\) 自身も有限。
各 image は三次元で、五 section の一次独立性から hyperplane に入らない。
\[
 5=L_t^3=\deg(f_t)\deg(f_t(Y_t))
\]
と五の素数性から、image は quintic hypersurface、\(\deg(f_t)=1\) である。
従って全 \(f_t\) は finite birational、smooth \(Y_t\) はその正規化となる。

これを scheme として family に戻す。
total image \(D\subset\mathbb P^4\times T\) は integral horizontal prime divisor。
ambient は regular なので Cartier であり、\(T\) を affine open of \(\mathbb A^1\) と取れば
\(D\) は homogeneous degree 五の一本の equation \(F_t\) で書ける。
horizontal prime equation に \(t-t_0\) 因子は無いので、各 \(F_{t_0}\) は零多項式ではない。
これは実際の irreducible quintic image \(f_{t_0}(Y_{t_0})\) を含み、同じ degree 五であるため、
**scheme として**その image に等しい。
特に \(D/T\) は flat な hypersurface family で、\(\omega_{D/T}\simeq O_D\)。
flat な正規化又は conductor の fiberwise 振舞いを先験的に仮定したのではない。

## 4. relative conductor と特殊化を同時に扱う

有限 birational \(f:Y\to D\) に対し、局所的に \(A=O_D\)、\(B=f_*O_Y\) と書く。
\(D\) は Gorenstein、\(Y\) は smooth over \(T\) で Cohen–Macaulay。
finite equidimensional map の下で \(B\) は maximal Cohen–Macaulay \(A\)-module であり、
\[
 \operatorname{Ext}^{i}_A(B,A)=0\quad(i>0),\qquad
 \operatorname{Hom}_A(B,A)\simeq f_*\omega_{Y/T}.                  \tag{3}
\]
後者は有限双対性と \(\omega_{D/T}=O_D\) である。
前者は Gorenstein 局所環上の MCM duality であり、\(B\) が locally free であるとは仮定しない。
MCM 性は、\(A\) の parameter ideal の \(B\) への拡張が parameter ideal となり、
CM な各有限局所 \(B\)-algebra 上で regular sequence となることでも確認出来る。

共通 fraction field で
\[
 \operatorname{Hom}_A(B,A)
   =\mathfrak c=\{b\in B:bB\subset A\}\subset A\subset B,
\]
同定は \(\varphi\mapsto\varphi(1)\)。従って (3) は自然な conductor inclusion
\[
 \omega_{Y/T}\longrightarrow O_Y,
 \qquad s\in H^0(Y,\omega_{Y/T}^{-1})                            \tag{4}
\]
を与える。

この map は**全 fiber の conductor inclusion と一致する**。
実際、\(\tau=t-t_0\) に対し \(0\to A\xrightarrow{\tau}A\to A/\tau\to0\) を dualize すると、
(3) より
\[
 \operatorname{Hom}_A(B,A)/\tau
  \simeq\operatorname{Hom}_{A/\tau}(B/\tau B,A/\tau).              \tag{5}
\]
evaluation at one もこの同型と可換であり、\(Y_{t_0}\to D_{t_0}\) は上で既に finite birational と証明した。
ゆえに (4) の特殊化は正確な conductor map で、中心で消える余分な \(t\) 倍を持たない。
必要なら \(B/A\) の \(T\)-flatness も、全 fiber の \(A/\tau\to B/\tau\) が injective なことから従う。
normalization-flat 又は conductor-basechange を未証明の入力にしていない。

## 5. 全 conductor-divisible tensor の descent と ambient lift

非零 \(\alpha_B\in H^0(B,\Omega^1_B)\) を固定し、(4) と contraction して
\[
 P_Y=i_{p^*\alpha_B}s\in H^0(Y,\Lambda^2T_{Y/T})                 \tag{6}
\]
を得る。中心では (1) の非零線を張る。
従って \(\alpha_B\) を一つの定数で調整すれば、\((P_Y)_0\) は**指定された** \(\pi_{Y_0}\) に等しい。
各 fiber で、三次元 bivector の Frobenius 条件はその対応一形式の integrability。
局所的に (6) の対応一形式は scalar 倍の closed \(p^*\alpha_B\) なので、全 \(P_{Y_t}\) は Poisson である。

\(\mathcal V=\Lambda^2T_{\mathbb P^4}\) とする。
\(df(P_Y)\) は \(f^*\mathcal V\) の section で、局所係数は全て conductor ideal \(\mathfrak c\) に入る。
\(\mathfrak c B\subset A\) より、これは \(\mathcal V|_D\) の正則 section \(\xi\) に descent する。
smooth \(D_t\) 上では \(f_t\) は同型なので \(\xi_t\) は tangent Poisson tensor である。
この descent は scalar のみ、又は normalization 上の rational tensor だけを lift したものではない。

quintic \(D/T\) の exact sequence は
\[
 0\longrightarrow\mathcal V(-5)=\Omega^2_{\mathbb P^4}
  \longrightarrow\mathcal V\longrightarrow\mathcal V|_D\longrightarrow0.
\]
\(H^0(\Omega^2)=H^1(\Omega^2)=0\) だから、relative direct image の同じ消滅により
\(\xi\) は \(\mathbb P^4\times T\) の regular relative bivector \(\Pi\) へ**一意に** lift する。
その中心は元の \(\pi_0\) である。
\([\Pi,\Pi]\) は全 smooth \(D_t\) 上で零であり、dense open から integral \(D\) 上でも零。
しかし
\[
 H^0(\mathbb P^4,\Lambda^3T(-5))=H^0(\mathbb P^4,\Omega^1)=0
\]
なので \([\Pi,\Pi]=0\)。つまり構成したのは**元の ambient Jacobi 方程式を全 parameter で満たす**実曲線である。

## 6. 実 Pfaffian と閉成分への全点包含

\(\Pi_t\) は \(D_t\) に tangent なので、その actual Pfaffian は \(F_t\) を割り切る。
両方 degree 五より
\[
 \operatorname{Pf}_{\rm actual}(\Pi_t)=c(t)F_t.
\]
中心は指定の rank 四 tensor だから \(c(0)\ne0\)。
零を残して \(T\) を再び縮めれば \(c(t)\) は unit、全 fiber で \(D_t\) は実際の被約既約退化 quintic である。
この操作は元の source の追加条件を課さず、各中心を失わない。

\(t\ne0\) では \(E_t\) は stable rank 三 degree 五。
受理済み genus 一認識 §§6／別作者review §6 により、任意のこの polarized normalization は
既知 \(q_{5,2}\) の normalization と同型であり、完全五 section の image は T40 の退化因子である。
さらに (1) と ambient restriction injectivity により、(6) の tensor 自身も projectively T40 である。
従って全 \(t\ne0\) の \([\Pi_t]\) は閉成分 T40 に入り、閉性により中心 \([\pi_0]\) も T40 に入る。

この証明は任意の \(B,A,F\) と任意の完全 section frame に適用する。
全 bundle moduli に一つの universal family を置く仮定は不要で、**各任意点**を同じ閉成分に結ぶので全 source の包含となる。
bundle による不安定型以外の追加 localization、normal R 仮定、clean 仮定、選択点や一般 source への制限は無い。
零 normalized tensor 又は actual Pfaffian 零は一般 rank 四という本定理の対象外であり、
それらをこの曲線から全 rank 二 support の新認識として主張しない。

## 7. 保存範囲

今回の新結論は **不安定 elliptic 全 rank 四 source の T40 包含**。
既に canonical／零 summand／cusp が空、stable は T40 と証明済みなので、
degree 五 rank 零・normalization sectional genus 一分岐の全 rank 四点が T40 に入る。
これは旧 G の像24以下を全 E と取り違えたのではなく、全 section／conductor／tensor を lift する新しい矢印である。
他の degree、rank 二 restriction、他の actual factor 型、全 E 又は総 component 数は主張しない。
新 component の増分は零。独立幾何review前の producer proof として保存し、旧受理済み証明は不変に保つ。
