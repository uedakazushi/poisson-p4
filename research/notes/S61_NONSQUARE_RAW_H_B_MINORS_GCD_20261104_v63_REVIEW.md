# B7×3の35最大minorとgcd1：独立監査

対象は `RAW_H_B_MINORS_COST`、`RAW_H_B_35_MINORS`、`RAW_H_B_MINOR_GCD` の凍結三packetである。sourceは既に独立受理されたraw H checkpointの同じ908項・7×3 blockであり、新しいGaussianやB-block探索は行わない。

## 元21成分と既知単位前処理

行と列のdivisorを `a,d,F,B,C,J` の非負整数冪と非零QQ定数から独立に作り直し、保存されたpolynomialと照合する。三つの非零QQ column multiplierも保持する。全21成分で

```
row_unit_i * prepared_ij * column_unit_j
  = QQ_multiplier_j * actual_source_ij
```

を独自のFraction疎多項式積で確認する。これらは `FULLU=D(adFBCJ)` 上の単位だけであり、pや未知のdeterminant因子を単位へ追加しない。prepared入力は整数345項である。

## 全35恒等式と元の二determinant

lexicographicな35三行組について、compact整数行列のdeterminantを独自の3×3 Laplace積和で求める。保存された全35値、元source minorへのrow/column scale換算、named factorとprimitive polynomialの全掛け戻しを確認する。QQ scalar倍だけによるdedupも独立に再現する。

零は7個、非零は28個、prepared minorの総項数は3525である。localized同じmaxminor idealのmonic generatorは22個、総2963項、最大次数15となる。さらに最後の二minor（35組中のindex33,34）は、compact値からscale換算するだけでなく、元908項source blockから直接determinantを再計算する。対応するprimitive polynomialは88項・次数10、および199項・次数14である。

このsource-boundな二多項式に対し、QQ上のexact multivariate gcdを独立に計算し、1を得る。したがって全22 generatorのgcdも1である。著者のnative executableやworker、GB、finite-field評価、点、factorization、primary decompositionは実行しない。

既存の407項Deltaも元の選択順で検証する。そのrow順は `(0,5,2)`、column順は `(0,2,1)` であり、両方の置換符号を保持する。review checkerの初期のmetadata assertionがcolumn順を誤って標準順に限定していたため、算術の一致を確認した後にそこで停止した。元の順序を正しく用いるよう修正し、行列式そのものの不一致とは区別した。停止logと修正前checkerはignored workに保持する。

## 正確な幾何学的結論

QQ[a,e,d]とその `adFBCJ` による局所化はUFDである。gcdが1の二generatorを含むidealを含むheight-one primeは存在しない。従ってBのrankbad locusは、FULLU上で空、またはcodimension≥2、すなわちdimension≤1である。unit ideal、空性、全FULLUでrank3はこれだけでは従わない。

この結論はB-blockについてのものに限る。既知のpairを用いるH5との還元はD(p)上に限られ、p=0で他の行がrankを回復するかは別課題である。全H7生成、full cubic57、quartic239、尽くしには昇格しない。別に保存された52秒GB timeoutも変更・再試行しない。

全入力hash、実測時間・メモリとexactチェック結果は専用review reportへ保存する。
