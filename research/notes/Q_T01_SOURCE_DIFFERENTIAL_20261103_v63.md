# Qa 全曲線上の actual T01 の特異性と source fiber

2026-11-03。対象は既約被約な **actual T01 成分**であり、Poisson
スキーム全体の接空間や、その全局所枝の分解とは区別する。

## 1. 結論

物理的な同次二次 bivector

\[
 Q_a=a x_0^2\partial_1\wedge\partial_4
 -\frac29x_0^2\partial_2\wedge\partial_3
 +\frac43x_0x_1\partial_2\wedge\partial_4
 +\left(\frac23x_0x_2+x_1^2\right)\partial_3\wedge\partial_4
\]

について、次を得た。

1. 全ての有限複素数 \(a\) に対し、\(Q_a\) は actual T01 の特異点である。
   より強く、その **affine cone T01 の Zariski 接空間次元は30以上**である。
   projective T01 では29以上であり、成分自体の次元はそれぞれ26、25である。
2. 下記の行列表示における31次元 source からの微分は、\(\mathbb Q(s)\) 上で
   階数24である。二つの表示の微分像は各24次元、共通部分18次元、和30次元である。
   共通部分はちょうど strict \(GL_5\) 軌道の17方向と \(Q_a\) 曲線方向の和である。
3. 各表示における **実際の source fiber** は一般点で正則7次元である。
   その接空間は、5方向の polynomial generator change と3方向の stabilizer
   conjugation で張られ、重複は1方向である。

1 は閉集合への拡張を含む証明済みの結果である。2 と3の階数は有理函数体上の厳密計算に
基づく。source fiber の主張は一般点であり、全ての特殊 \(a,s\) に拡張した主張ではない。
**二つの \(s\) の値から、T01 の局所枝が二本あることや normalization fiber が二点ある
ことは導かない。** また微分階数24だけを理由に target の特異性を推測していない。

## 2. 実際の T01 への多項式写像

既に独立検証済みの
[`Q_T01_MEMBERSHIP_20261103_v63.md`](Q_T01_MEMBERSHIP_20261103_v63.md)
の写像を用いる。\(M\in\operatorname{Mat}_5\)、\(1\leq k\leq4\) に対して

\[
 W_k(M)=M^k-\frac{\operatorname{tr}(M^k)}5I,\qquad
 \Psi(M,c)=\sum_{1\leq i<j\leq4}c_{ij}(W_i(M)x)\wedge(W_j(M)x).
\]

全ての \(W_k\) は可換かつ trace-free なので canonical divergence と Jacobi が消える。
regular semisimple \(M\) では4本が Cartan の基底をなし、6個の \(c_{ij}\) は全ての
toric period を動かす。従って既約 source \(\mathbb A^{25}\times\mathbb A^6\) の像は
全体として actual T01 に入り、閉包はその affine cone である。ここでは singular
\(\beta\) に \(\Delta^{-1}\) を含む公式を代入していない。

既知の actual T01 の次元26は、regular semisimple 点における20次元の Cartan conjugacy
class と6個の period から得られる。projective dimension は25である。入力 catalogue
`legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json` の T01 項とも一致する。
source 全体が T01 に写ることは、今回の接空間を **reduced T01 自体**の接空間に
含めるための必要な依存結果である。

\(s\neq0,-3\) として、strictly lower triangular matrix \(N_s\) の非零成分を

\[
 (N_s)_{10}=1,\quad (N_s)_{21}=s,\quad
 (N_s)_{32}=\frac{s+3}{2},\quad (N_s)_{43}=-\frac9{s+3}
\]

と置く。\(N_s^5=0\)、\((N_s^4)_{40}=-9s/2\neq0\) であり、regular nilpotent である。
以下を取る。

\[
 a(s)=\frac{4s+6}{3s(s+3)},\qquad
 c_{14}=\lambda=-\frac{2a(s)}{9s},\quad
 c_{23}=\mu=-\frac4{9s^2(s+3)},\quad c_{ij}=0\ \text{otherwise}.
\]

このとき **全150係数で** \(\Psi(N_s,c)=Q_{a(s)}\) である。二次 cover

\[
 3as^2+(9a-4)s-6=0
\]

の他方の根は

\[
 \sigma(s)=-\frac{3(s+3)}{2s+3},\qquad a(\sigma(s))=a(s)
\]

である。以下の \(s\) と \(\sigma(s)\) の比較は有理函数体 \(K=\mathbb Q(s)\) 上で
行う。\(2s+3=0\) 等の特殊値にこの表示の式をそのまま代入する必要はない。

## 3. full150 微分と source fiber の正則性

source の列順は \(M_{00},M_{01},\ldots,M_{44}\) の25列、続いて
\(c_{12},c_{13},c_{14},c_{23},c_{24},c_{34}\) の6列である。行順は
\(0\leq i<j\leq4\) の辞書順、その内部で \(0\leq p\leq q\leq4\) の辞書順の
\(x_px_q\partial_i\wedge\partial_j\) である。

各 \(E=E_{hk}\) に対して Fréchet 公式

\[
 dW_j(E)=\sum_{r=0}^{j-1}N^rEN^{j-1-r}
 -\frac15\operatorname{tr}\left(\sum_{r=0}^{j-1}N^rEN^{j-1-r}\right)I
\]

を用い、積の微分から \(J_s=d\Psi_{(N_s,c)}\) の全 \(150\times31\) 成分を作った。
厳密 RREF の階数は24であり、従って fiber の embedding dimension は7である。

この7次元を実際に積分する。まず strict stabilizer

\[
 H=\{g\in GL_5:gQ_a(g^{-1}x)g^t=Q_a(x)\}
\]

を用いる。直接物理係数を微分した Lie action

\[
 L_GQ=GQ+QG^t-(Gx)\cdot\partial Q
\]

の全 \(150\times25\) 行列は階数17である。従って \(\dim H=8\)。標数0の代数群なので
\(H\) は滑らかである。計算ではさらに

\[
 L_GQ=J_s([G,N_s],0)
\]

を全150行、全25列で確認した。regular nilpotent の centralizer
\(C_{GL_5}(N_s)\) は5次元であり、\(H\) に含まれる。従って
\(g\mapsto(gN_sg^{-1},c)\) は実際の fiber 内に3次元の軌道を与える。

次に \(M'=p_0I+p_1M+\cdots+p_4M^4\) を、\((p_0,p_1,p_2,p_3,p_4)=(0,1,0,0,0)\)
の近傍で動かす。\(M=N_s\) では \(W_1(M'),\ldots,W_4(M')\) は
\(W_1(N_s),\ldots,W_4(N_s)\) から可逆な4行4列変換で得られる。この wedge-square の
逆変換で \(c'\) を選べば、\(\Psi(M',c')=Q_a\) が **正確に** 成り立つ。これは
5次元の実際の fiber 内変形であり、単なる first-order solution ではない。

コードは \(\dot M=N_s^r\), \(0\leq r\leq4\) と、それぞれ必要な \(\dot c\) を明示的に
作り、\(J_s(\dot M,\dot c)=0\) を確認する。これら5方向と \(H\)-conjugation の3方向の
和は階数7である。重複の1方向は \(\dot M=N_s^4,\dot c=0\) である。実際
\(G=(-9s/2)E_{41}\) とすれば \([G,N_s]=N_s^4\) であり、
\(I+tG\) の conjugation はちょうど \(N_s\mapsto N_s+tN_s^4\) を与える。

よって滑らかな \(H\) と generator changes の積から実際の fiber への射の微分階数は7。
fiber の局所次元は7以上であり、embedding dimension 7と一致するため、この fiber の
局所環は正則7次元である。これは **source fiber の正則性**であり、T01 の正則性や
source が normalization であるという主張ではない。一般の semisimple fiber の5次元と
比べて、ここには実際に追加の fiber directions が存在する。

## 4. 二つの微分像の和は30次元

\(J_s\) と \(J_{\sigma(s)}\) を同じ150係数空間に並べる。全
\(150\times62\) 行列の厳密 RREF 階数は30である。report に記録した行・列で選ぶ
30次小行列式は

\[
 \frac{2^{36}(2s-3)(2s+3)^4(4s+9)
 (2s^2+6s+9)^9(2s^2+12s+9)}
 {3^{32}5^6s^{19}(s+3)^{18}}\neq0.
\]

この非零小行列式だけで接空間次元30以上を証明するには十分である。階数の上限30と
各像24次元も RREF で厳密に確認したため、共通部分は18次元である。

共通部分には strict \(GL_5\) 軌道17方向が含まれる。また
\((N_s,c(s))\) を \(s\) で微分すると

\[
 J_s\frac{d(N_s,c(s))}{ds}=a'(s)\frac{dQ_a}{da},\qquad
 a'(s)=-\frac{2(2s^2+6s+9)}{3s^2(s+3)^2}\neq0
\]

となる。コードはこの恒等式を全150行で確認する。軌道の17方向と
\(dQ_a/da\) の和の階数は18であり、他方の表示の微分像にも含まれる。従ってこれが
共通部分の全体である。

補助的な特殊点確認として、\(a=1\) では
\(s=(-5\pm\sqrt{97})/6\)、\(a=2/9\) では \(s=(3\pm3\sqrt5)/2\) を取り、対応する
厳密な代数体上でも \(24,24,30\) を確認した。この二点確認を一般点の代用には使わない。

## 5. 全ての有限 a への結論

\(\Psi\) は actual reduced T01 に写るので、上記の二つの微分像は共に
\(T_{Q_a}\mathrm{T01}\) に含まれる。従って一般の \(a\) において

\[
 \dim T_{Q_a}\mathrm{T01}\geq30>26=\dim\mathrm{T01}.
\]

有限型既約被約複素 variety の特異点集合は閉じている。その \(a\mapsto Q_a\) による
逆像は \(\mathbb A^1\) の閉集合であり一般点を含むので、全 \(\mathbb A^1\) である。
これは cover の分岐値 \(a=\pm4i/9\) と \(a=0\) も含む。

さらに接空間次元が30以上となる集合自体も閉じているので、30という下界も全ての有限
\(a\) に成立する。全ての \(Q_a\) は係数 \(-2/9\) を持つ非零 bivector である。
affine cone の非零点の近傍では scalar \(\mathbb G_m\) 因子を分離できるため、projective
T01 では接空間次元29以上、局所次元25となる。

ここで用いた閉性は、二つの \(s\) の表示を分岐値まで別々に延長する必要をなくす。
これは T01 の局所枝の本数、unibranch 性、normality、completed ideal、あるいは全
Poisson 成分の列挙を決定しない。

## 6. 再現と証拠

新規 producer は
`code/q_t01_source_differential_20261103_v63.py`。次の一回の実行で、約7秒、通常の
Python/SymPy のみを使用して再現できる。M2 や巨大 GB は用いない。

```sh
python code/q_t01_source_differential_20261103_v63.py \
  --output work/q-t01-source-differential-v63/replay
```

保存先は `reports/Q_T01_SOURCE_DIFFERENTIAL_20261103_v63/` である。
`report.json` に全ての階数、各非零小行列式の行列番号と値、特殊点の代数体計算を記録した。
`matrices.json.gz` に全150行の両 Jacobian、直接計算した Lie action、stabilizer basis、
5個の generator change、conjugation directions、実際の fiber を張る方向、および
共通18次元空間を疎係数形式で保存した。行・列番号は全て **0-based** である。

入力 hash は membership producer、その frozen report、T01 catalogue、および今回の
producer を固定する。membership 自体の独立検証は
`research/notes/Q_T01_MEMBERSHIP_20261103_v63_REVIEW.md` にあり、今回の微分計算とは
別の証拠である。
