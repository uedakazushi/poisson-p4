# CA原点を通る実polynomial self-pair：独立検証

`AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63` を受理する。これは同じactual CA tensorを表示する、異なる **二つのCA旗** のQQ[a]上の全次数の族である。有限jetから延長可能性を推測したものではない。

## 物理tensorの別構成

raw CAのHamiltonian式を直接canonical化し、座標交換を行って全150×20 fixed-flag基底Cを再構成した。保存されたcritical packetの基底と一致する。さらに証明書の150係数を引用せず、次のvector fieldから第一tensorを構成した。

\[
D_a=(x^2+12ahx/5)\partial_h-3ax^2\partial_x/5
 +(hx+2axz/5)\partial_z+(h^2+z^2-8axw/5)\partial_w,
\qquad P_a=\partial_y\wedge D_a.
\]

全150係数で \(P_a=P_0+a(-3C_1-C_9)\) に一致した。独立なphysical式での二つの係数minorはrows16,20で−12/5であり、非零な係数row20は全aで−1である。従ってa=0を通るordinaryな非定数projective tensor arcであり、途中で零tensorにならない。

## 全次数の二source membership

\(M=E_{10}-2E_{03}\) について、\(M^3=0\)、\(M^2\ne0\) を確認し、

\[
g=1+aM+a^2M^2,\quad g^{-1}=1-aM,\quad\det g=1
\]

をQQ[a]で直接掛け合わせた。次に、producerの疎多項式ringでの外積計算とは別に、5×5のphysical行列を用いて

\[
Q_a=g\,P_a(g^{-1}X)\,g^T
\]

を全次数で展開した。全150係数と全20 source係数が保存値に一致し、\(Q_a=C\alpha_2(a)\) を確認した。\(\alpha_2(a)\) はaの二次以下のQQ多項式であり、aや未知因子を逆転しない。逆向きのfinite pushforwardも直接計算し、\(g^{-1}_*Q_a=P_a\) に戻ることを確認した。

第一・第二tensorのdivergenceとstrict cubic Jacobiを、physical変数による直接のsymbolic微分で検証した。両方の25 divergence係数と350 Jacobi係数はQQ[a]で恒等的に零である。producerの新規再実行による証明書のbyte一致も確認した。

## 同じtarget tensorの二つの旗

固定旗を \(F_0=(L=\langle\partial_y\rangle,W=\langle h,x\rangle,H=\langle x\rangle)\) とする。\(P_a\) と \(g_*P_a\) はどちらもfixed-flag CA sourceに属するので、同じ \(P_a\) の二表示は \(F_0\) と \(g^{-1}F_0\) である。

covectorへの正しい逆変換を行い、後者を独立に計算した。

\[
L'=\langle\partial_y\rangle,\quad
W'=\langle h-2az,\ x+ah-2a^2z\rangle,\quad
H'=\langle x+ah-2a^2z\rangle.
\]

二つのW'基底のh,x minorは1、L'とのannihilationも零であり、全aで実際の旗になる。H'とxの外積のh,x成分はaである。従って **全a≠0で二旗は異なる**。単なるsource tangentやparameterizationの差ではない。

全族はconstant vector factorを持つためrank4ではない。Pfaffianは零である一方、\(E\wedge P_a\) の012成分は \(x^2(x+3ah)\) で全aに対して非零多項式なので、各aでgeneric projective Poisson rankは正確に2である。

## 範囲

このQQ[a]恒等式だけで、ordinary CA self-pairが任意次数で存在しないという主張は否定される。一方、全proper fiber数、二sourceのgeneric immersion、CA原点の全odd ideal、全conductor、全枝や接合型は本packetから分類しない。CA/CB間の異なる旗を転用した構成でもなく、Nの別source模型へ二つのCA旗を同一視したものでもない。

独立codeは `code/ambient_CA_origin_self_pair_review_20261103_v63.py`、receiptは `reports/AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63_REVIEW.json` である。旧四次jet witnessはsource発見の由来として保存されているが、今回の証明はそのjetの外挿に依存しない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/ambient_CA_origin_self_pair_review_20261103_v63.py \
  --replay work/CA-origin-polynomial-pair-producer-replay \
  --output work/CA-origin-polynomial-pair-independent-review.json
```
