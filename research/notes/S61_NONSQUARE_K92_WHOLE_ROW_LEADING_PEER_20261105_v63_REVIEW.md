# K92 row-cleanup の非零係数を直接五角 border で独立検証

2026-11-05。元の二対二 block による fused cleanup を再利用せず、保存済み ABC72/D49 の literal 係数から

    det [[A4, B_col0], [C_row2−59375 C_row0, D_20−59375 D_00]]

の一つの係数を検証した。各入力は同じ整数分母 L で払う。五行の subset determinant DP を独立に実装し、元の v,u,w prefix や主計算の積恒等式を使わない。

先に非負の絶対係数で support・高さを見積もった。要求次数6663、同時 DP support は最大204、係数335 bit、積769対以下、全 owned 上界63272、入出力 RAW 約3.02MB以下となった。60秒・2GB・150000 scalar と metadata・32MiB の制限内である。

符号付きの一回の実計算は695対、322 bitとなり、L五乗で戻した6663次係数は

    17149069629000000000 / 7 ≠ 0

で、主計算の同じ一係数と一致した。従って、この formal Kronecker specialization 上で全 row-cleanup の比例関係は成立しない。検証範囲はこの非零有理係数であり、主計算の全4131係数一致やそれより前の全零係数を独立に確認したとはしない。K92全十一角 minor の非零、全FULLU57、quartic239の結論ではない。

generic rational callback 自体を univariate 多項式で置換してはならない。generic への適用は別の finite-source と a-potential lift の証明を経由する。非零 specialization は generic callback の非零を示せるが、零 specialization から generic callback の零は従わない。

再現は同名 Python の `--mode source`、続いて `--mode value`。出力先を新しい空フォルダとし、既存の frozen receipt を上書きしない。
