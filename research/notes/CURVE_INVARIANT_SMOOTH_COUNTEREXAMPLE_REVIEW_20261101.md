# 滑らかな canonical curve 上の保存則反例：独立点検

2026-11-01。`CURVE_INVARIANT_SMOOTH_COUNTEREXAMPLE_20261101.md` の反例を承認する。
Δ≠0 であっても、全 S61 上で c4,c6,j が基点移動により一定になるとはいえない。

## 1. 大域的 Poisson 構造

六つの log-canonical 係数 (2,8,12,3,5,7) を持つアフィン bivector を、
多項式自己同型 `y0=z0+z1², yi=zi (i>0)` で押し出して独立に再構成した。
結果は原ノートの同次二次行列

\[
 B_{01}=-2x_2^2,\qquad
 (B_{12},B_{13},B_{14},B_{23},B_{24},B_{34})
 =(2x_1x_2,8x_1x_3,12x_1x_4,3x_2x_3,5x_2x_4,7x_3x_4)
\]

の射影アフィン化と一致する。
従って単なる局所的な非射影座標変換の例ではなく、全 P4 上の二次射影 bivector である。
押し出しの構成に加え、アフィン係数を直接微分して Jacobi 恒等式も確認した。
アフィン稠密開集合で Jacobi が零なので、大域的射影 Jacobi も零となる。

## 2. 独立な二点の正規化

基点 (t,1,1,1)、t=2,3 において、元実装の Darboux frame と primitive-jet 関数を使わず、
別の symplectic Gram–Schmidt を実装した。得た frame F について `F P(p) Fᵀ=J` を確認し、
一次微分を変換して trace を除いた。

この別経路の primitive β に、既に点検された Fisher 不変量 utility を適用すると次を得る。

| t | c4 | c6 | Δ | j |
|---:|---:|---:|---:|---:|
| 2 | 6371584 | 15372597248 | 12934592726630400 | 3946956213086464/197366222025 |
| 3 | 3347089 | 5894745263 | 6364426689994275/4 | 149989816182953863876/6364426689994275 |

原ノートの値と全て一致する。両方で Δ≠0、j は相異なる。
同一の大域的 P に沿う不変量は有理関数なので、異なる二値が得られることだけでも
非定数性を証明する。標数零であるためその微分も恒等的には零でない。

従って提案された全方向の保存則は Δ≠0 上の全 S61 では成立せず、対応する輸送条件の
全てを `radical(I_S61:Δ^∞)` に属するものとして追加することはできない。
この反例だけでは、さらに小さな一般一次部分の開集合上の可能性は決まらない。

## 3. 点検範囲

新しい実装は `code/curve_invariant_smooth_counterexample_review_20261101.py`、
証明書は `reports/CURVE_INVARIANT_SMOOTH_COUNTEREXAMPLE_REVIEW_20261101.json`。
押し出し、射影への大域化、Jacobi、二点の frame・primitive β を別途再構成した。
Fisher tensor 縮約の実装自体は既存の検証済み utility を共有しており、今回別実装にしていない。
QQ(t) 上の c4,c6 の巨大な普遍的縮約も再実行していない。
T02 への装飾図式の同定と像次元15は、このレビューの対象外である。

```sh
.venv/bin/python code/curve_invariant_smooth_counterexample_review_20261101.py \
  --output reports/new-run/smooth-counterexample-review.json
```
