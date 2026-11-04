# B18 chain lift：最初の60秒試行と境界表現の診断

2026-11-04、v63。最初の typed operator 試行は60秒で停止した。
source とサイズの検査までを保存し、\(UP=QB\) はこの試行から主張しない。
最大観測 RSS は約764 MBであり、新しい physical source solve や高次 jet を行ったものではない。

停止後の read-only 依存調査で、既知 scalar 関数に対する負号境界の欠落を特定した。
例えば CA graph64 の既知 odd 関数は node17870=`n(13367)` である。
元 DAG はその負号を簡約して直接 node13367 を参照する。
最初の reader は17870だけを既知関数として止め、13367を元 source の中へ展開した。
このため不要な symbolic 展開が増えた。既存 HONEST10 で既に必要だった負号閉包と同じ問題である。

最初のコード・hash・失敗 report・最後の progress は同名 report directory に保存した。
別途承認された改訂試行は
`AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63`
であり、既知関数 \(i=-j\) の境界に \(j=-i\) も加える。
これは関数を変更する操作ではない。最初のコードや証拠を改訂版で上書きしていない。
