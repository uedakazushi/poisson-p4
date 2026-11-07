# rank零plane・m二・rankM二・M(p)非零：全像dimension二十三以下

2026-11-07、主義務G。producer exact arithmetic及び全source幾何接続。
独立review前は新global claimへ混ぜない。

**命題。** `PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_FIRST_JET.md` の全八familyは
projective Poisson像dimension≤23を満たす。
M(p)=0、rankM三、O(0)及びm≥3はそれぞれ別の射程であり、本結果だけで全plane又は全Gは閉じない。

## 1. 全sourceと最弱linear upper bound

前掲noteは全Lax systemを、NDSS(s≠0)、NDSSN、NDNSS、NDN、
BSS(τ,s), Δ=4τ−3s≠0、BSSR(τ)、BNS(τ)、BNR(τ)へ尽くす。
全normal/base/tensor変換とその逆を明示し、BSSのresonant rowを別のBSSRへ保持した。
rankM二かつM(p)非零という全supportを、normal-form anchorの標本で置き換えていない。

全30 normal-weight一Π1、全9 normal-weight二Π2を保持し、
全350 projective Schouten係数からD1及びD2を作る。
全higher Poisson fiberは

`d≤dimkerD1+dimkerD2`

を満たす。fixed Π1でΠ2の方程式はD2のinhomogeneous systemであり、
残る非線形higher Jacobiはこのsupersourceを狭めるだけだからである。

## 2. reusable runnerと有限契約

`code/classification_plane_rank_zero_fullrank_weights.py --rank2-nonkernel` を使った。
全旧FunctionDef ASTは、immutable rankM二/M(p)零reportに埋め込まれたproducer bytesと
runtimeで同一と確認した。旧source・report・review bytesを変更していない。
jobは `jobs/plane_rank_zero_rank_two_nonkernel_weights.json`、gate PASS後に120秒/1536MiBで実行した。

証明書は `reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_NONKERNEL_WEIGHTS_20261107.json`。
全八familyの350×30、350×9、actual six-shear30×6、完全kernel、exact rank minorを保存する。
実行は30.939秒、最大RSS70,248KiB。

初回は契約通りone-variable factor `2s²−1` で停止した。
そのimmutable runtime reportが当時のproducer bytesを保持する。
有限algebraic rootを捨てず、同一source・同一operatorへexact number-field rank演算だけを追加し、
新bindingで再gateして上記全計算を終えた。
二parameterのnonlinear curveを無証明でgeneric扱いするfallbackは使っていない。

## 3. 一parameter及び固定sourceの全結果

| whole source | generic dimkerD1 | generic dimkerD2 | 全generic image upper | 全finite wall image upper |
|---|---:|---:|---:|---:|
| NDSS(s), s≠0 | 6 | 0 | 23 | 22 |
| NDSSN | 6 | 0 | 22 | — |
| NDNSS | 6 | 0 | 22 | — |
| NDN | 7 | 0 | 22 | — |
| BSSR(τ) | 6 | 0 | 22 | 22 |
| BNS(τ) | 6 | 0 | 22 | 21 |
| BNR(τ) | 7 | 0 | 22 | 23 |

NDSSのselected-minor productは非零constantを除いて

`s⁹(s−1)(s+1)(3s−1)(3s+1)(2s²−1)(3s²−1)(18s²+1)(72s²+1)`。

全rational root及び四つのirreducible quadraticの全rootを保持した。
quadratic factor fについてはexact field QQ[α]/fで全matrixのrank/minor/kernelを検査した。
保存kernelへの全行代入はこのfieldで厳密に零、minorは非零。
QQ上のirreducibilityから全complex conjugate embeddingでrankは等しいので、
一つのfaithful field embeddingがfの全rootを同時に証明する。
数値的な√近似又は一方のrootだけの推測は用いない。
s=0は算術上保持したがm≥3なので明示的に除外する。

BSSRの全minor wallは `τ=−3/2,1/6,3/2,1/2,−1/2`。
BNSは `τ=1,−1/2,−3/2,1/6,3/2,1/2`。
BNRは `τ=1,−3/2,−1/2,3/2,1/2`。
BNRの1,3/2,1/2ではd=9、−1/2ではd=8、他のnonexcluded wallではd≤7である。
fixed τによりbase dimensionが一つ下がるので、最大image boundは二十三である。

## 4. BSSの全二parameter curve・point cover

formal reportではetaが幾何noteのsを表す。
QQ(τ,s)上のgeneric kernelはD1が六、D2が零。
selected-maximal-minor productは全十七linear factorに分かれた。

`2τ+1, 2τ+3, 3s+2, 2τ−3, 2τ−1, 3s−2τ+3, 3s+4τ−4,`
`6s−10τ+3, 3s−6τ−1, 3s−4τ−2, 3s−4τ+1, 3s−2τ−3,`
`3s+2τ−1, 6s−6τ+1, 3s−2τ−1, 3s−2τ+1, 3s−4τ`。

最後のΔ=0 lineはrow零のm≥3補集合であり、本BSS sourceから除外する。
その上のm二row非零は、上のBSSR whole familyとして別に保持済み。
残る十六lineの各々をexact rational parameterで全てparametrizeし、
その**全350行**について新しいQQ(r) maximal minorとcomplete kernelを計算した。
全line genericでd=6である。
各line自身のminorの全rootをさらにexactに検査し、合計142個のpoint certificateを保存した。
このうち14個はΔ=0へ入る明示補集合、残る128個は全てd≤8である。

line相互のintersectionをgenericとして落としていない。
各lineの全pointは、そのline minorが非零のrank-certified patch又は同minorの全rootに入る。
したがって全intersectionもこの二者のどちらかで覆われる。
global generic matrixの一つのminorを検査するだけで止めたものではない。

危険なd=8のnonexcluded pointは、重複を除いて

`(τ,s)=(3/2,4/3), (1,5/3)`

だけである。
他のnonexcluded line-rootはd≤7。
generic source dimension十二、lineで高々十一、pointで高々十なので、
plane六・projectivization一を戻して

`generic: 12+6+6−1=23`,
`line: 11+6+6−1=22`,
`point: 10+8+6−1=23`

以下となる。全二parameter BSSも上界二十三で閉じる。

## 5. 全base parameterの復元と結論

全八sourceのsymbol base dimensions12/11/10と、fixed oneparameter rootでの一dimension減少は
前掲geometry noteで導出した。
本reportのoperator rankはgradingを保つ実可逆normalizationで元の全sourceへ戻る。
全表と二parameterのcurve/point coverを合わせ、全projective像は二十三以下。
よって未知dimension≥25 componentの一般元は、このwhole first-jet stratumに入らない。

K又はnormal coneのsymbol因子だけから、全higher Dが可約だとは推論しない。
全Jacobi nonlinear support、個々のknown component帰属又は全E点分類も主張しない。
O(1)のrankM三・rankM二/M(p)零はそれぞれ別証明で≤23/22。
rankM一はglobal-vector-field型のsupersourceへ入れられ、O(0)型とともに次の未完了課題である。
plane multiplicity三・四及びrank二restrictionの残余も保持する。
