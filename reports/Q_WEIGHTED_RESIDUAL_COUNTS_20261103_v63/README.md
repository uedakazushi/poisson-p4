# GF(32003) 上の固定 a 全残差展開と有界な次元探索

`a=2/9`（有限体では 7112）、正の重みの自由変数34個について、保存済み35変数有限モデルの全350残差を native FLINT で展開した。jet に切り詰めていない。

- 全350式中228式が非零、総18,041,430項。
- 全74 implicit pivot 式は多項式として零。重み1–5の残差はすべて零。
- 全350展開を独立した有限体 DAG 評価と照合した。
- 展開124.168秒、最大RSS 2,341,400 KB、最大中間多項式1,864,073項。参照数による解放、6,000,000,000 byte / 300秒の上限を使用した。
- 重み8以下の22式（26,291項）の完全GBは150秒で打ち切った。次数12までの部分GB34要素の先頭単項式イデアルの次元は29。包含 J ⊂ in(I_weight≤8) ⊂ in(I_full) により、固定 a の全有限体モデルの次元は29以下とだけ結論する。
- 全350式の四種類の Lex 先頭単項式イデアルの次元は32,32,33,32で、上界29を改善しない。

全モデルの正確な次元・Hilbert級数・既約成分分解は未決であり、この記録では標数零への次元移行を主張しない。35変数の完全な有限モデル自体の証明とは区別する。

`report.json`、`residuals.summary`、`weight-counts.json` が全展開の小さな保存記録である。469MBのM2入力と1.3GBの疎係数表は ignored work 内にあり、下記ソースと公開済みモデルから再生成できる。`initial-bound-*` と `lex-*` は有界計算の実行ソース・結果を保存している。巨大な部分基底520001項の文字列保存は行わず、先頭単項式だけを記録した。

## 再現

未使用の work 出力名を選び、リポジトリの既設 Cloud 環境で実行する。

```sh
source /workspace/cloud-setup/activate.sh
python code/q_weighted_residual_expansion_20261103_v63.py --a 2/9 --output work/q-full-fresh
python code/q_weighted_residual_hilbert_20261103_v63.py --expansion work/q-full-fresh --output work/q-w8-fresh --weight-cut 8 --prepare-only
python code/q_weighted_residual_initial_bound_20261103_v63.py --partial-run work/q-w8-fresh --output work/q-bound-fresh --degree-limit 12 --timeout 90
python code/q_weighted_residual_leading_scan_20261103_v63.py --expansion work/q-full-fresh --output work/q-lex-fresh --timeout 90
```

次数12の再現計算は32秒で終了した。部分基底の保存処理で終了した旧probeは ignored work に残し、採用結果は小さなsummaryへの直接出力による最終成功版とした。

## 独立した固定 a 四変数確認

`fixed-a-four-coordinate-validation.json` と `fixed-a-four-coordinate-polynomials.json` は、a=2/9で weight6 の z11,z12,z23,z33 だけを残した FULL350 閉部分sliceの独立QQ検証記録である。u=z11−4t, v=z12−6t, w=z23+4t, t=z33 とすると残差イデアルは (u,v,w)(u,v,w,t)。9×9係数行列式534600、全350 physical Jacobiとの一致を確認した。

これは既存の [a を自由にした強い閉部分slice定理](../../research/notes/Q_WEIGHT6_CLOSED_SUBSLICE_20261103_v63.md) の固定 a における独立検証として保存し、新しい別定理とは扱わない。再現コードは `code/q_weighted_resonance_four_slice_20261103_v63.py`。
