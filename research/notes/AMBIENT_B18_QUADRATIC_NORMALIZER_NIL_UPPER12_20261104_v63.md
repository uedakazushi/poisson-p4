# B18 実交叉の nil 指数の無条件上界12

2026-11-04。対象は実交叉 (D=CA\cap CB\) の既定 B18 点 (P_0\) における
形式的完成環である。全 Poisson スキームや全大域交叉への主張ではない。
二つの実 quadratic normalizer と既受理の whole (E_4\) を用い、
新しい source・jet・点・CAS 計算を行わずに

\[
  \operatorname{Nil}(\widehat{\mathcal O}_{D,P_0})^{12}=0,
  \qquad \operatorname{Nil}(\widehat{\mathcal O}_{D,P_0})^2\ne0
\]

を示す。従ってこの局所環の nil 指数は3以上12以下である。三乗が零かどうか、
正確な指数および有限長 nil-cube の長さは決めない。原点の深さ19の新しい結論から
nil 三乗の非零性を推論することもない。

## 1. 二次拡大の任意の基底変換核

補題。(A\subset B\) は可換環の包含で、(B=A+A\tau\) および

\[
 \tau^2=s\tau+\nu,\qquad s,\nu\in A
\]

を満たすとする。2が (A\) で可逆なら、任意の可換 (A\)-代数 (M\) に対し

\[
 Q=\ker(M\longrightarrow M\otimes_A B),\qquad Q^2=0.
\]

証明。自由多項式表示 (M=A[X]/J\) を取る。変数集合が有限である必要はない。
核の元 (k\in A[X]\) は (JB[X]\) に属する。

\[
  k=j_0+j_1\tau,\qquad j_0,j_1\in J
\]

と書けるのは、任意の (B\) 係数を (A+A\tau\) に分解できるためである。

ここで (A[X]\to B[X]\) は包含である。従って (a=k-j_0\in A[X]\) について

\[
 a^2-sj_1a-\nu j_1^2=0
\]

という (B[X]\) 内の等式は (A[X]\) 内にも成り立つ。右二項が (J\) に入るので

\[
 k^2\equiv a^2\equiv0\pmod J.
\]

最後に (2kk'=(k+k')^2-k^2-k'^2\) により核 ideal 全体の平方が零となる。
基底変換後の包含、flatness、自由な (A\)-basis (1,\tau\) は仮定していない。

特に trace-free な (s=0\) では (a^2=\nu j_1^2\) だけで足りる。
この証明は nil 元の個別の平方零性だけで止めず、核 ideal の平方零性を与える。

## 2. 完成基底変換と実二段比較

実零次数 source は共通 physical regular 底 (R_0\) 上で

\[
 A_i=R_0\oplus I_i\tau_i\ \subset\ B_i=R_0\oplus R_0\tau_i,
 \qquad \tau_i^2=\nu_i\in R_0,\quad i=A,B.
\]

この同じ (A_i,B_i\) と同じ zero ambient (P_0\) が、受理済み
`ZERO_NORMALIZER_TENSOR_MODEL` の source と独立監査により指定されている。
二族の trace basis は同一視しない。(R_0\subset A_i\) より (B_i=A_i+A_i\tau_i\)
である。各拡大は有限だが (A_i\)-flat とはしていない。

\[
 D_0=A_A\widehat\otimes_{P_0}A_B
   \longrightarrow D_1=B_A\widehat\otimes_{P_0}A_B
   \longrightarrow E_4=B_A\widehat\otimes_{P_0}B_B
\]

とする。第一写像は (A_A\subset B_A\) の (D_0\) への基底変換、第二は

 (A_B\subset B_B\) の (D_1\) への基底変換である。有限 module の completed tensor
は通常の tensor の完成であり、ここでは結果が有限 complete module なので、
補題を各段階にそのまま適用できる。別法として有限 topologically free power-series
表示で証明しても、(J\) が閉じているため同じ係数分解と恒等式が成立する。

各段階の核を (Q_1,Q_2\)、合成の核を (K\) とすれば

\[
 Q_1^2=Q_2^2=0,\qquad K^2\subset Q_1,\qquad K^4=0.
\]

合成の像を第一環と同一視したり、tensor の Tor を捨てたりしていない。
第二環で (K\) の像の平方が零なので (K^2\subset Q_1\) となるだけである。

## 3. 正規化比較から上界12

受理済み proper-source CI と五支持枝の独立監査は、embedded localization も含む
**whole (E_4\)** について

\[
 \operatorname{Nil}(E_4)^3=0
\]

を証明している。被約商や core のみを代用していない。
環写像は nilpotent を nilpotent に送るので、(N_0=\operatorname{Nil}(D_0)\) に対し

\[
 N_0^3\subset K,\qquad N_0^{12}=(N_0^3)^4\subset K^4=0.
\]

この上界には generic conductor の支持分類、(K\) の rank、深さ、(T\) の長さ、
または比較写像の Spec-surjectivity は不要である。実二次包含と whole (E_4\)
の nil-cube 零性だけを使う。

## 4. 正次数、群因子、下界

既受理の相対非負 graded algebra を先に有限 (w\)-polynomial 環として局所化すると、
actual stable tail (M\) を得る。(w=0\) 完成環に (w=1\) を代入しない。
全正次数の係数写像 (D_n\to M\) は (T_1=T_2=T_3=0\) と安定化により単射であり、
whole tail の nil ideal は平方零である。従って nil homogeneous 元の積は、
総次数が正ならすべて零となる。零次数と正次数との混合積も含む。

nilradical は閉じた torus-stable ideal なので、有限 Vandermonde 射影と完備性により
各 homogeneous 成分も nilradical に属する。積の恒等式は有限商を通して完成に移り、
任意の (q\ge2\) について

\[
 \operatorname{Nil}(D_{\mathrm{trans}})^q=N_0^q.
\]

従って全横断完成環にも上界12が成り立つ。共通 actual (G_{19}\) の regular
power-series 因子、および射影 chart の scalar 因子を戻し／除いても同じ nil 指数である。
power-series 環の nil ideal は元の有限生成 nil ideal の延長であり、その被約商は
被約環上の power-series 環となるためである。

一方 actual (Q_b\) family は (P_0\) を通り、その非空 generic open の whole D 環で
nil 二乗が非零である。nil-square ideal は coherent なので支持は閉であり、
その閉包上の (P_0\) でも非零となる。忠実 flat 完成と実群 product はこれを保存する。
これは既受理 nil-cube finite-origin corollary とその独立監査が明示した下界である。
従って

\[
 3\ \le\ \text{nil index at the B18 germ}\ \le\ 12.
\]

上界8を得るには、さらに **全域** の (N_0^2\subset K\) を証明すればよい。
幾つかの generic source 模型で像の平方が零であるだけでは、その全域包含は従わない。
本稿は E4 の全 nil 支持との source 同定や hidden finite-origin 部分を消したとはしない。
上界8、三乗零、三乗非零および正確な指数は未主張である。

依存 path・SHA256、補題の普遍多項式恒等式と指数の bookkeeping は同名 report に保存する。
