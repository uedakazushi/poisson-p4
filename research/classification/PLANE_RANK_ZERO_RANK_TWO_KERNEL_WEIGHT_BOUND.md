# rank零plane・m二・rankM二・M(p)=0：全像dimension二十二以下

2026-11-07、主義務G。producer exact arithmeticと全source幾何接続。
独立review前はglobal完了主張へ混ぜない。

**命題。** `PLANE_RANK_ZERO_RANK_TWO_KERNEL_FIRST_JET.md` の全stratumは
projective Poisson像dimension≤22を満たす。
coefficient map rank二だがM(p)≠0、rank三、O(0)、m≥3は本命題の外に保持する。

## 1. 一parameterを含む全四source

前掲noteはactual base/normal変換とtensor scalar及びその逆を明記し、
全sourceをND、BSSρ（全ρ∈C）、BN1、BN0へ尽くした。
BSSρを手で0又は1へ置換せず、nilpotent rowの零/非零だけを実変換で二orbitにした。
全normal weight一30及びweight二9は、rank三の場合と同じcanonical Euler quotientで保持する。

全350 projective Schouten coefficientで
`D1=E∧[Π0,Π1]`、`D2=E∧[Π0,Π2]` を作る。
actual higher Poisson fiberは、Π1の全D1 kernelと、その各点でのD2 inhomogeneous systemに入るので

`d≤dimkerD1+dimkerD2`。

高weightの非線形Jacobiはこれを狭めるだけである。

## 2. 旧証明bytesを保持したreusable runner

`code/classification_plane_rank_zero_fullrank_weights.py --rank2-kernel`
というisolated modeを、rank三の独立review freeze後に追加した。
元の全FunctionDef ASTは、immutable rank三reportに埋め込まれたproducer bytesと
runtimeで完全一致することをassertした。
旧report、旧review、旧proofのbytesを更新していない。

新jobは `jobs/plane_rank_zero_rank_two_kernel_weights.json`。
分類gate及びinput hashesを確認してから実行し、120秒/1536MiB以内、
4.744秒・最大RSS56,808KiBで終了した。
証明書は `reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_KERNEL_WEIGHTS_20261107.json`。
全350×30、350×9、actual six-shear30×6、complete kernel、nonzero rank minorを保存する。
BSSのformal parameter名tauは本noteのρを表し、reportにもこの同一視を明記する。

## 3. 全kernelと全一parameter例外

| whole source | dimkerD1 | dimkerD2 | d upper |
|---|---:|---:|---:|
| ND | 6 | 0 | 6 |
| BSS一般ρ | 6 | 0 | 6 |
| BSSρ=3,1,−1 | 6 | 0 | 6 |
| BSSρ=−3 | 7 | 0 | 7 |
| BN1 | 6 | 0 | 6 |
| BN0 | 7 | 0 | 7 |

BSSのselected-maximal-minor productは非零constantを除いて

`(ρ−3)(ρ−1)²(ρ+1)²(ρ+3)`。

全四rootをQQ上で評価し、全行でkernel/minorを検査した。
ρ=0はgeneric patchに含まれ、未処理の零parameterではない。
generic-rank claimだけでexceptional sourceを捨てていない。
いずれのD2も全sourceでinjectiveである。

## 4. actual full sourceのdimension復元

前掲noteのfixed-plane affine symbol dimensionsはND十一、BSS一般十一、
BSS固定ρ十、BN1十、BN0九である。
plane choice六を足し、tensor projectivization一を引くと

| whole source | 全projective image upper |
|---|---:|
| ND | 11+6+6−1=22 |
| BSS一般ρ | 11+6+6−1=22 |
| BSS固定ρ=3,1,−1 | 10+6+6−1=21 |
| BSS固定ρ=−3 | 10+7+6−1=22 |
| BN1 | 10+6+6−1=21 |
| BN0 | 9+7+6−1=21 |

となる。
gradingを保つ実normalizationにより、このrank結果は各whole orbit/familyへ戻る。
fixed normal formのkernel dimensionを、消費したbase parameterを復元せず使ったものではない。
generic openとfinite四rootでBSS全ρを尽くし、全四sourceで上界二十二を得る。
従って未知dimension≥25 componentの一般元はこのstratumに入らない。

## 5. symbolのactualDと残余

全weight零symbolのactualDは

`ND: 4bc*(bz²−cw²)`、
`BSSρ: −6bcw*(2bz+cw)`、
`BN0,BN1: 2b²w*(2bz+cw)`。

全sourceに非線形higher coefficientsを加える前の因子型であり、
これだけから全higher sourceのD因子を推論していない。
本命題は全30/9必要operatorによる弱いdimension boundで完結する。
all nonlinear Jacobi、個別成分帰属又は全E点分類は求めない。
他のplane first-jet strata及び全Gは別の未完了義務として残る。
