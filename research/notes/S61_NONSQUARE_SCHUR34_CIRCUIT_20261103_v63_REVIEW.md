# 第34段階の完全 Schur 回路の独立検証

## 結論と有限データ

`S61_BINARY_NONSQUARE_SCHUR34_CIRCUIT_20261103_v62` の全回路を独立に受理した。新しい行列は厳密な `23×58` の多項式行列であり、その全1,334成分を有限の算術回路で与えている。展開済み native cache は第33段階の `24×58` のままであるが、回路としては第34の構成的な単位 pivot が完了している。

独立コード `code/s61_nonsquare_schur34_circuit_review_20261103_v63.py` は producer の `build` / `evaluate` を import しない。保存した各 node を底から形式的な式の署名に翻訳し、全6,244 node の種類・arity・参照先・位相順序、全1,392個の元行列参照、33個の非零 cofactor 参照、24個の完全な image 及び全1,334根を検証する。全 node が出力から到達可能であり、余分な未検証演算はない。元行列の410MBの hash はストリームで再確認する。

元行列 `A` と58成分の `q` は `S61_NONSQUARE_AFTER33_SCALAR_UNIT_20261103_v63_REVIEW` の独立整数演算で受理済みのものに hash で一致する。特に `A_7q=g=aP²J` は探索や有限体の予測ではなく、QQ 多項式環の恒等式である。native の第3変数は `P=ad` を意味する。

残す元の行は `[0,...,6,8,...,23]`、列は `[0,...,57]` であり、先頭7行が積の行、残る16行が genuine cubic の行である。各 node の検証により厳密に

\[
v_i=\sum_{k=0}^{57}A_{ik}q_k,\qquad
B_{ij}=gA_{ij}-v_iA_{7j}
\]

を得る。零の `q_k` の項だけを省略しており、次数による切捨てや jet への置換はない。

## Cokernel の証明

任意の可換環で `g` が単位であるとする。`v=Aq` を列として追加しても列 image は変わらない。pivot 行7を保ち、他の行を `gR_i-v_iR_7` に変える行変換の行列式は `g²³` で単位である。その後、各旧列から `(A_7j/g)` 倍の追加列を引く列変換は行列式1を持つ。行と列を並べ替えると正確に

\[
\begin{pmatrix}g&0\\0&B\end{pmatrix}
\]

となる。従って `coker(A)≅coker(B)` である。反転するのは既知の `g=aP²J` だけで、cofactor や未知の core 因子は反転しない。この議論は非被約の係数環でも成立する。

## 別の多項式 fixture

作者の例とは異なる QQ[s,t,e,h] 上の例を使い、`s` だけを反転する。

```
A = [[1,0,0],[0,0,h],[st,-s,0]], q=(e,te-1,0), pivot row=2.
```

ここでも非零の cofactor はどちらも指定した単位ではない。`A_2q=s` であり、独立計算した結果は

```
B = [[s(1-te),se,0],[0,0,sh]].
```

行・列変換そのものを構築し、行列式 `s²` と1、正確な単位 block を検証した。元の行列式は `sh`。結果の最大 minor は `0,s²h(1-te),s²he` であり、第2 minor と `t` 倍の第3 minor の和が `s²h` となる。従って minor ideal は厳密に `(s²h)` であり、未知因子 `h` は消えていない。

## 再現と範囲

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/s61_nonsquare_schur34_circuit_review_20261103_v63.py \
  --output work/schur34-independent-review
```

保存した report は `reports/S61_NONSQUARE_SCHUR34_CIRCUIT_20261103_v63_REVIEW/report.json`。入力 hash、全参照の件数と別 fixture の結果を含む。大きな外積の展開や巨大な新 cache の作成は必要ない。

構成的な厳密回路の pivot 数34と、展開済み native cache の pivot 数33を区別する。`23×57` への列置換、全点での product41 / full cubic57 / quartic の被覆、ambient exhaustion はここでは主張しない。
