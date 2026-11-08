# reduced R・m三の全九source：独立Fraction算術監査

2026-11-08、`/root/quartic_sources`。producerとは別作者のexact replay。
全九必要fullhigher sourceの算術及びdimension/factor含意を受理した。
root-coverの全幾何的尽くしは別作者の別reviewへ分担し、本稿の算術から推論しない。

## 1. 再構成した全データ

対象は `reports/CLASSIFICATION_PLANE_REDUCED_M3_NINE_FULLHIGHER_20261108.json`。
embedded producer sourceのSHAを照合したが、producer又はそのarithmetic moduleをimport/実行していない。
独立 `classification_plane_m4_special_review.py` のFraction/sparse5 engineを再利用した。
CAS、浮動小数、有限体、numerical rank又は補間は使っていない。

二first-symbol JSONの全raw bracketsと全parallel target basisを独立AST parserで読み、
各profileで次を元定義から再構成し、保存matrixの全entryと比較した。

* 全150 P0、元全350 projective Jacobi、実Euler-Pfaffian零。
* 全350x30 J1、350x9 J2。
* 実D3の全24x30 map（base degree二、normal degree三の全24monomials）。
* 全leading-target coefficient rows、complete annihilator、annihilator-times-D3を追加した全必要J1。
* 全150x39 Euler/Levi/scale group matrix。

全kernel certificateは、nonzero maximal minor、rank+nullity、全free-coordinate identity、
元の全rowへのkernel代入零で完全性を確認した。
target annihilatorは独立RREFからのcomplete kernelとも一致した。
全30/9 basisは独立engine自身がraw36/Euler6 complementを確認する既source定義を用いる。

## 2. 全sourceへのdimension/factor含意

全actual normal symbolが指定parallel targetへ入るなら、P1は全target-constrained ker J1へ入る。
固定P1に対してP2はweight二Jacobiの空又はker J2 affine-torsor部分集合なので、
全higher fiber dimensionはk1+k2以下。
未使用のweight三/四Jacobiは解を減らす方向にしか働かない。
leading amplitude零を許した大きな必要sourceでもこの結論は有効。

Euler columnsのrankは24。
全13Levi directions及びP0 scaleを加えたmatrixのrankを24+Bとすると、
projective Levi orbit dimensionはB-1、Levi stabilizerは13-B。
全六filtered base-from-normal shearsも同じassociated-graded P0上のhigher sourceを保つため、
fullPGL parameterizationのfiber dimensionは19-B以上となる。
従って全projective像は

\[
 24-(19-B)+k_1+k_2=B+k_1+k_2+5
\]

以下。各Bは入力のfirst-symbol上界とも一致した。
whole parallel-target spaceはP0のintrinsic connectionから定まりstabilizerに保存される、という幾何的入力をここで使う。
任意のchosen target sliceにこのgroupを流用する意味ではない。
その全parallel/root-cover scopeは別幾何reviewと合わせて読む。
dimensionを本質的に使うM/SF_C3_UNUSEDのtargetは、その全intrinsic parallel spacesである。
他の七profilesの結論は全u因子なので、そこでgroup上界を使わずにも空性が従う。

factor型ではP0の全u-rowがu倍、全constrained K1のbad v^2 u-rowが零であることを独立確認した。
P2にはnormal columnsが無いので全solutionのu-rowがu倍。
u零のprojective tensorは三次元hyperplaneにtangentとなりrank二以下、従って実Dはuで割れる。
これは全必要sourceへのproper factorで、normal leading symbolだけの因子を全actualDへ外挿していない。

## 3. 全九結果

| profile | B | constrained k1 | k2 | 粗像上界 | 指定sourceの結論 |
|---|---:|---:|---:|---:|---|
| S | 7 | 8 | 3 | 23 | actualu因子、既約Dで空 |
| N | 9 | 8 | 3 | 25 | actualu因子、既約Dで空 |
| D | 9 | 9 | 3 | 26 | actualu因子、既約Dで空 |
| M | 10 | 7 | 0 | 22 | 全像22以下 |
| SF_C3_UNUSED | 9 | 8 | 1 | 23 | 全像23以下 |
| TRIPLE_HX2 | 9 | 8 | 3 | 25 | actualu因子、既約Dで空 |
| TRIPLE_HY2 | 9 | 6 | 1 | 21 | actualu因子、既約Dで空 |
| TRIPLE_HYZ | 9 | 9 | 3 | 26 | actualu因子、既約Dで空 |
| TRIPLE_HXY | 9 | 7 | 1 | 22 | actualu因子、既約Dで空 |

25/26という粗数値を未閉componentの存在と呼ばず、全factorで指定locusを空にする。
M又はSF_C3_UNUSEDの個別component帰属/存在は、この上界だけから加えない。

## 4. 実行receiptと旧bytes

新isolated `--reduced-m3` modeを既reviewer
`code/classification_plane_quintic_casimir_review.py` に追記した。
旧10 FunctionDef ASTは全て不変。
直前exact reviewer及び全二current契約bytesは
`inputs/plane_reduced_m3_nine_review_frozen_inputs_20261108.json` に保存した。
旧reportは改変していない。

独立review契約は実行前gate PASS。
一回の全九source再構成は1.049秒、RSS25,708KiBでALL PASS。
reportは `reports/CLASSIFICATION_PLANE_REDUCED_M3_NINE_SEPARATE_REVIEW_20261108.json`。
ここで証明したscopeは全九必要sourceの算術とそのG非優越/空性である。
別作者の `PLANE_REDUCED_NONCONCURRENT_M3_ROOT_COVERS_SEPARATE_REVIEW.md`
（SHA 77103082f373d9b25b39a9dc82878cdb16920184e9f0f285d3bb8b3357ba2b70）は
横断空性及び接線四/triangle五の全root-cover必要尽くしを独立に受理している。
この別scope receiptと接続すれば、対応する全nonconcurrent reducedR/m3 G branchへ用いられる。
concurrent及びnonreducedRを含む全m3又は全分類の完了へ拡張しない。
