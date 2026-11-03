# 新しい一般係数模型の有界実行と、残る計算制約

2026-11-03。既存セットアップは再実行せず、導入済みCloud M2 1.24.11と
既存Singular 4.4.1 runtimeを再利用した。
ここに記す実行は、以前の18・9・8パラメータGBの無変更再実行ではない。
新しい87次定数拡大、および六パラメータmonic判定を具体的な厳密入力にした。
**全S61のunit、一般アフィン有限性、五点尽くしは得ていない。**

| 模型 | 係数・未知数 | 制限 | 結果 |
|---|---|---|---|
| 87次定数拡大の通常多項式表示 | QQ(s5,c1,c2)[z,q45]、z^87−2、元134式と二線形式 | 180秒、仮想メモリ4 GiB | timeout 124、GB未完 |
| monic六変数関数体 | QQ(s5,c)[q44]、一線形式でq32を正確に消去、元134式 | 180秒、仮想メモリ4 GiB | memory cap 14、GB未完 |
| cを通常多項式変数にするblock消去 | QQ(s5)[q45,c]、元134式とmonic式、q45を先に消去 | 120秒、仮想メモリ3 GiB | timeout 124、GB未完 |

メモリ値は設定したvirtual memory上限で、実RSSのpeakと混同しない。
Singularは全て `--cpus=0` で、並列modStd workerは起動していない。
実行は `timeout --signal=TERM --kill-after=10s` と `ulimit -v` で区切った。
M2は `M2 --script INPUT`、Singularは
`work/nonsquare-v62/singular-runtime/run-singular.sh --cpus=0 -q INPUT` を使った。

通常多項式版を試した理由は、六変数関数体版が初期係数消去で4 GiBへ達したためである。
cを係数体から外し、元134式のコンパクトな形を保って、
qについてのblock orderで非零のcのみのeliminantを求める別の計算にした。
非零eliminantが得られればk(c)上のunitと同値であり、元generatorとのcofactor積を確認する。
この検査自身は二座標toyで正確に通過し、eliminantはc^4−6c^3+9c^2−2だった。
初回の未完GBは、後処理の未対応 `deg(poly,poly)` を短い別試験で発見したため中断し、
`diff(poly,var(i))` に修正してtoyを確認した後、表の最終120秒実行を行った。
最終実行で後処理には到達していない。

別にQQ(α)(s,c)というnested係数体でのM2 GBはengine非対応だった。
通常多項式表示はその問題を回避し、同じ関数体上の小さなGBと実change matrixは成功した。
従って全実行の未完を、環境全体の未導入やセットアップ失敗とはしていない。

producer `code/s61_affine_new_projection_attempts_20261103_v62.py` は、
各実行の厳密入力・出力・終了receiptをgzipで保存し、
凍結した二つの数学的判定producerとの入力一致と、第三模型の全134行保存を確認する。
数学的なunit結論を表す `all_passed` は使用せず、receiptの検証成功を区別する。

次に必要なのは、真に小さい係数模型、厳密な構造的消去、または別の正当化された計算法である。
同じ模型をそのまま延長する場合は、上の明示的な時間・メモリ上限を超える計算資源が必要となる。
既存の数値点の空性を一般unitへ移すことや、未確認のp整有限性を仮定することはしない。
