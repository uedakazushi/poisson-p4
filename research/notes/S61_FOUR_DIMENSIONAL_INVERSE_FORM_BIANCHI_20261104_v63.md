# 4次元 inverse-form 恒等式は Pfaffian 倍 Bianchi と Koszul の和

2026-11-04。残る43方向の候補として、非退化 bivector の逆2形式から
syzygy を作る方法を検査する。結論は、以下の操作自体は新しい primitive
source を与えず、**既存 Bianchi の Pfaffian 倍に明示的な Koszul 関係を
加えたもの**になる、という普遍的恒等式である。一般の非単項式飽和、
actual135 の全代入、四次 rank239 は未解決のまま保つ。

## 1. 規約と多項式恒等式

標数0の微分多項式環で、4座標を `x0,...,x3` とする。

```
       [ 0  a  b  c ]          [ 0 −f  e −d ]
P  =   [−a  0  d  e ],   A =  [ f  0 −c  b ],   p=af−be+cd.
       [−b −d  0  f ]          [−e  c  0 −a ]
       [−c −e −f  0 ]          [ d −b  a  0 ]
```

ここで `PA=AP=p Id`。`p!=0` 上では、行列 `A/p` を係数とする2形式を
`omega` とする。2形式の係数は `i<j` に対して記し、同じ係数を全添字で
書く場合には通常の `1/2` を付ける。

```
J_ijk = sum_l(P_il partial_l P_jk + P_jl partial_l P_ki
                                      + P_kl partial_l P_ij),
j_i = (−1)^i J_(0,...,hat(i),...,3),
X_j = sum_i partial_i P_ij.
```

3形式にも同じ dual 規約を使う。このとき、**分母を含まない**恒等式は

```
j = −grad(p) − A X,
E = p dA − dp wedge A = p² d omega,
dual(E) = p X + P grad(p) = −P j.                 (1)
```

`p² d omega` という表示だけが `D(p)` 上の記法であり、`E` 自体と
以下の恒等式は `p=0` を含む多項式恒等式である。さらに

```
dE = 2 dp wedge dA,
p dE − 2 dp wedge E = 0.                         (2)
```

従って逆形式を微分して得る式は、`d²=0` から自動的に成り立つ。
問題は、それを Jacobi source の cofactor vector として書いたときの内容である。

## 2. Bianchi との差を cofactor で計算する

実際の `j` と区別して、形式的な独立 vector `h` を導入する。
凍結済み full Schouten operator と同じ符号では、4-vector の係数は

```
B_P(h) = X dot h − sum_ij P_ij partial_i h_j.     (3)
```

実際の Jacobi dual を入れると `B_P(j)=0` が Bianchi 恒等式である。
一方、(1) の `−Pj` の場所に `−Ph` を入れて (2) を書いた operator は

```
D_P(h) = p div(−Ph) + 2 grad(p)^T P h.
```

Leibniz 則だけで

```
D_P(h) − p B_P(h)
  = 2(grad(p)^T P − p X^T) h
  = −2 j^T P h.                                 (4)
```

最後の等号は (1) と `AP=p Id`、`A^T=−A` を使う。
`h=j` とすると `j^T P j=0` だから、(2) は **`p B_P(j)=0`** に戻る。

より強く、(4) の差は cofactor vector の段階で

```
−2 sum_(i<j) P_ij (j_i h_j − j_j h_i)
```

という明示的な Koszul 関係である。単に二つの多項式が零になることを
比べたのではない。actual Jacobi source を `j=Q F` と書けば、差の cofactor は
`−2 F^T Q^T P Q` となり、`Q^T P Q` の交代性により元の `F` の標準
Koszul 関係の多項式結合になる。この主張は `Q` が正方・可逆でなくても成り立つ。

従って空間単項式ごとに係数を取っても、この inverse-form construction は
既存 Bianchi の係数関係と Koszul 関係の加群内に入る。Pfaffian の空間係数は
通常の convolution multiplier にすぎない。actual135との source correspondence
自体は既存の full operator 監査を参照し、今回はその巨大な再展開を行わない。

## 3. Pfaffian で割る段階に残る正確な gate

(4) の差を除く前の cofactor が `p` で割れるわけではない。例えば block-diagonal
`P01=a,P23=f` で、`−2 j^T P h` の `h0` 係数は `−2 a² partial_1 f`。
これは普遍多項式環で `p=af` の倍ではない。したがって、式が `D(p)` 上で零と
いうだけで cofactor vector の分母を消したことにはならない。

Koszul を引いてから割れば、(4) は単に `p B_P(h)` なので、得るものは元の
Bianchi に戻る。この特定の操作から新しい43本を得たとはいえない。

一般の飽和を排除する結論でもない。有限な actual source coefficient ring `S` で、
既知 syzygy vector の加群を `M subset S^N`、局所化する既知因子を `g` と書けば、

```
M[g^−1] intersect S^N = M : g^infinity.
```

新しい候補 `v` を逆形式・分母付き計算から採用するには、少なくとも
`g^k v in M` と `v in S^N` を **全 source cofactor について**確認する必要がある。
新方向を主張するなら、さらに `v notin M` を示す必要がある。
`M:g^infinity=M` は今回証明していない。localization 後の等式、個々の source
多項式の消失、ある点での値はその代用にならない。

ここで空間変数に依存する `p(x)` の逆元を導入する場合は、まず
`S[x]` の cofactor module でこの contraction を行い、最後に空間係数へ戻す。
`p(x)` を根拠なく `S` 内の一つの単元と見なさない。
既存の一単項式除算 preflight と parameter-field constant-combination negative は
限定された class の結果であり、この一般 contraction gate を解決しない。

## 4. 一回の小さい普遍検証

専用コード:
`code/s61_four_dimensional_inverse_form_bianchi_20261104_v63.py`。
報告・全式:
`reports/S61_FOUR_DIMENSIONAL_INVERSE_FORM_BIANCHI_20261104_v63/`。

6個の0-jet、24個の1-jet、60個の対称2-jet、形式 `h` とその微分20個からなる
110変数の `QQ` 多項式環で、(1)–(4)、(2)、actual Bianchi、Schouten の符号、
Koszul の交代零を全係数で検証した。9項目すべて PASS、実測0.086951秒、
最大RSS48,076KiB。60秒・2GB上限内の一回だけの universal symbolic test である。

新しい beta 点、有限体、GB、full135 expansion、quartic matrix、rank search はない。
結果は普遍的な polynomial/operator identity の証明と、その小さい全係数検証であり、
full80 kernel・残る43方向・full239・一般 beta exhaustion の解決ではない。

再現時は未使用の `--out` を指定する。

```
.venv/bin/python code/s61_four_dimensional_inverse_form_bianchi_20261104_v63.py \
  --out /tmp/inverse-form-bianchi-replay
```

参照した凍結済み frontier は報告の `input_sha256` に記録した。
従来の本文・報告・source bytes は変更していない。
