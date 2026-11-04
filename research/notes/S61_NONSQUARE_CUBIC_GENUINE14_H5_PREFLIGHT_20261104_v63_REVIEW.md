# Genuine14/H5のliteral bindingと既知単位因子：独立監査

対象は `S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63`。前段で独立受理されたfinal14×21の同じcokernel表示と、raw H-linearの固定checkpointを結ぶ。著者producerは実行せず、新しいGaussian、単位探索、B-blockのminorやGB、標本、rank probeを行わない。

## 70根を元checkpointへ戻す

残る真正source labelは `[10,12,13,14,15,16,17,18,19,20,21,22,23,24]`、heavy target labelは `[27,28,29,30,31]` である。final14×21の各根を、`matrix_parent` から親の `checkpoint` へ直接追い、元の行・列を全70箇所で一致させる。零根も元の零多項式へ対応することを確認する。この二段のliteral pathにはPOWER ancestorが無く、次数207のpure係数は一つも展開しない。

## 単位因子と全70恒等式

係数環は `QQ[a,e,d,1/(adFBCJ)]`。F,B,C,Jは既存sourceの同じ多項式であり、未知のpやQは単位に含めない。固定transcriptのinitial divisorと一回の完了pivotだけからτ_iを読み、`rho_i=lambda_i*tau_i/s`、`s=F B² C/12500000` を単位の指数と有理定数で復元する。逆に `lambda_i=rho_i*s/tau_i` の分子・分母を作り、元のsource clearing polynomialそのものへ掛け戻して照合する。因数探索や新しい除算scanではない。

各row unitは非零有理数倍のC、aFC²、またはFC²であり、column unitは最後のpsi列だけ `gamma=-4FBC²/625` である。これらはすべてFULLU上の単位である。全70の旧checkpoint分母は定数なので、独自の疎QQ[a,e,d]多項式演算で

```
old_numerator = old_denominator * rho_i * compact_ij * gamma_j
```

をそのまま検算する。元70成分の13854項を保持し、compact側2407項・非零43・最大148項と各列の415/304/302/302/1084項を確認する。heavy部分が零のsourceは正確に10,13,16であり、そのpure tailを零とは呼ばない。

L-coordinate correctionの相殺も普遍的な可換多項式恒等式として確認する。従ってこのheavy blockの同一視は正当なfull-source row rescalingとpsi target rescalingへ延長でき、未展開のpure tailを保持する。raw15×6の既受理のGaussian transcriptを再実行する必要はない。

## 階数の範囲

前段の独立raw-H generic-rank gateと全70単位同値により、actual heavy14×5は一般係数体でrank5である。既受理の明示chartは `FULLU∩D(pQ)` に限る。p=0やQ=0でこの一つのminorが零でも、actual rankが落ちるとは結論しない。全FULLUでのheavy rank5、H7生成、full cubic57、quartic239、尽くしや新しいkernel数は本監査の結論に含めない。

実行時間・メモリ、全依存hashと恒等式結果は専用review reportに保存する。失敗して保存された旧representation比較のsource/log/receiptは変更しない。
