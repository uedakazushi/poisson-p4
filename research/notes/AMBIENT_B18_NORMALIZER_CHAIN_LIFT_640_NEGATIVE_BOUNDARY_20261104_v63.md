# B18 実 normalizer への640列 chain lift

2026-11-04、v63。全て同じ actual 共通 (R_0) 上での有限 operator 表示である。

既受理の実行列

\[
 P:R_0^{640}\longrightarrow R_0^9,\quad
 U:R_0^9\longrightarrow R_0^4,\quad
 Q:R_0^{48}\longrightarrow R_0^4
\]

に対し、具体的な \(B:R_0^{640}\to R_0^{48}\) を構成し、

\[
 UP=QB
\]

を全640列で確認した。(B\) の非零 scalar operator 成分は2808個、
typed scalar DAG は17,849 nodes、証拠は914,583 bytesである。
実 scalar 係数の巨大な多項式展開、元の (r/g\) の再除法、追加 source solve や jet は行わない。

## 1. 実 source と二つの C0 基底を保持する

source は既受理の `FULL_RELATIVE_INPUT`、`T0_HONEST_C0`、
`ZERO_NORMALIZER_TENSOR_MODEL` と同じファイル・node ID である。
全 input hash と両 source prefix の一致を保存した。
scalar 関数参照は C0 source の**node ID**を用い、同じ表示名の変数を同一視しない。
古い正則 (r/g\) 関数は元の関数を表す opaque 参照としてそのまま再利用する。

実 honest C0 基底は

\[
 (1,z_{A,0}\tau_A,\ldots,z_{A,6}\tau_A,g_A\tau_A).
\]

一方、保存済み612列の Q0 の乗数は

\[
 (1,z_{A,0}\tau_A,\ldots,z_{A,6}\tau_A,f_{A,63}\tau_A).
\]

最後の乗数を (g_A\tau_A\) へ置き換えない。
両者の正確な関係
\(f_{A,63}=g_A+\sum e_{A,i}z_{A,i}\) は既受理の actual CA prefix-Hadamard 証拠を使用する。
今回はこの関係を、実関数 node と係数 node の全てを保持したまま利用した。

## 2. 十二 physical target 上の68個の CB 多項式

weight zero の CB 方程式68個は、28個の線形 Koszul、36個の二次積、四つの線形 graph である。
target は実 physical rows

\[
 18,19,28,29,60,67,72,63,64,68,73,74
\]

の十二座標であり、source の正則 scalar 係数を伴う target 次数二以下の多項式である。
CB の最初の八 odd 変数を、既存の unipotent 行列により
\(Z_0,\ldots,Z_7\) へ替え、係数を \(b_0,\ldots,b_7\) とすると

\[
 b_iZ_j-b_jZ_i,\qquad Z_iZ_j-\nu_Bb_ib_j
\]

が最初の64個である。ここで (b_7=g_B\)、
\(Z_7=Y_7-\sum e_{B,i}Y_i\) という**元の actual source 変換**を保持する。
四 graph の係数も既存 CB regular membership の元関数を使う。

各多項式 (f_k\) に対して

\[
 f_k(q_A)-f_k(q_B)=\sum_{j=1}^{12}(q_{A,j}-q_{B,j})H_{k,j}
\]

という具体的な target-Hadamard を構成した。
線形の場合はその係数、積の場合は

\[
 H_{ij,t}=L_{i,t}Z_j(q_A)+Z_i(q_B)L_{j,t}
\]

を用いる。ここで (L\) は上の actual unipotent target 変換の勾配である。
全68個について、保存済み actual CA 評価対と再構成した評価対が一致し、
上の Hadamard 恒等式が有理数上の形式多項式として成り立つことを確認した。
CB 側の消滅には、既受理の正確な (g_B\) と四 graph の membership 恒等式を使用する。
閉点で零になる関数を単位として反転しない。

## 3. 四成分表示と chain identity

正確な normalizer tensor は基底
\((1,\tau_A,\tau_B,\tau_A\tau_B)\)、
\(\tau_A^2=\nu_A,\tau_B^2=\nu_B\) の (R_0\)-自由四加群である。
各 (H_{k,j}\) をその四成分で保存し、上の**古い actual C0 乗数九個**を掛ける。
得られた四成分を、実 (Q\) の column 順序 (4j+c\) に並べたものが (B\) である。

最初の28列は C0 の intrinsic Koszul なので零 lift とした。
残り612列の各々について、次を確認した。

1. 元 CB 多項式と乗数の順序が、保存済み actual Q0 の source 対に一致する。
2. exact honest C0 の係数表示を (U\) で送ると、その実 source 対になる。
3. 実 (Q\) の四十八列をそのまま用いた (QB\) が同じ四成分になる。

従って、同値な別 presentation を仮定せず、保存済み行列に対して実際に (UP=QB\) が成り立つ。
source の有限関数恒等式を代入した後の scalar 式は、元の正則関数を自由係数として扱う
可換有理多項式の正規形で検査した。

## 4. 負号境界の正確な扱い

既知 scalar 関数 (i=\mathrm{n}(j)\) の境界では、子 (j\) にも \(-i\) を割り当てる。
CA/CB の q0/full0 odd 関数十六組について、元 source の有理線形算術でその別表示が同一であることも検査した。
全負号境界の node ID・子 node ID・対応 operator を、計算の最初に別ファイルへ保存した。
これにより、DAG の \(-(-b)\) 簡約から元 source の内部へ不用意に展開しない。
未束縛の normalizer source 変数へ降りた場合は即時停止する。

最初の60秒試行はこの境界処理が欠けて停止し、別 packet に保存してある。
今回は同じ上限60秒・2 GBで承認された別表現による試行であり、1.888秒・約215 MBで終了した。

## 5. 完了した範囲

`actual-UP-QB-640-typed-operator.json` は全68 Hadamard、全640 source/lift 対応、
(48\times640\) の sparse (B\)、既存 (U,Q\) の四成分、全関数参照を含む。
既受理の complete normalizer syzygy operator により、これを regular8＋零4 の source へ
制限すれば次の低次数比較 complex の ancestry を組める。
その unit 行列の制限・比較 homotopy・特殊化 kernel の計算は、この packet 自体では行っていない。
比較 kernel (K\)、元の深さ、一般 DVR rank に関する新しい結論はない。

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_normalizer_chain_lift_640_negative_boundary_20261104_v63.py
```
