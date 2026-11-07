# rank零plane・globalVの二つのK零残余を一つのJacobi条件で閉じる

2026-11-07、主義務G。exact characteristic-zero producer証明。
対象は`PLANE_RANK_ZERO_GLOBAL_VECTOR_LOW_PARAMETER_WEIGHT_BOUND.md` §3で
linear upper二十五を残した二つのwhole source：

- normalss、T=J,V=J²,K=0。
- normalnil、T=J+J²,V=J,K=0。

いずれもfixed-plane affine first-symbol base dimensionは九であり、
full weight一D1のkernel八、weight二D2のkernel三であった。
Π1を全八係数u0,…,u7で書き、全Π2九係数を保持する。

## 1. full weight二の必要条件

`Q2(u)=(1/2) E∧[Π1(u),Π1(u)]`を全三百五十係数vectorとする。
weight二Jacobiは厳密に

`D2*r+Q2(u)=0`。

reportのD2 rank六のnonzero rational minorを使い、full D2を消すcokernel functional φを構成する。
`φD2=0`を全九columnで検証し、`φQ2`が非零polynomialであることを証明した。
したがって任意のfull Poisson solutionで`φQ2(u)=0`。
これは特定のhigher pointにおける不成立から全sourceの不存在を推論する議論ではなく、
全source上で必須の一つの非零polynomialを得る議論である。

uはreportに保存したfull30×8 kernel basisの係数で、三百五十rowの番号はzero-basedである。
normalssについて

`φ=(2/3)e51−e61+e132`,

`φQ2=−4*(u2−u7)^2/9`。

normalnilについて

`φ=(1/3)e51−(2/9)e52+(2/9)e53−(1/3)e61+(1/9)e62−(2/27)e63+e131`,

`φQ2=−2*(9u1²−15u1u3+36u1u5−12u1u7+5u3²−30u3u5+5u3u7+36u5²−24u5u7−u7²)/9`。

両者とも明らかに非零であり、C上の全zero locusはdimension高々七。
weight二equationの各非空fiberはkernel D2のtranslateなのでdimension三。
従って全higher solutionはdimension高々十となる。
weight三・四Jacobiはこのsourceをさらに制限するだけである。

## 2. ambient像と証明書

actual全normalizationで得たsymbol source九、plane choice六、tensor projectivization一を戻すと

`dim image≤9+7+3+6−1=24`。

したがって両whole K零sourceはGで未知優越成分の一般元になれない。
K非零sourceは前のfull operator packetで既に二十一以下なので、二つの一parameter族全体が閉じる。
全known-component membership又はその閉包が既約成分であるとは主張しない。

immutable証明書は`CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_ONECUT_20261107.json`。
full eight-parameter kernel、全Q2三百五十係数、D2 minor inverse、φ、非零polynomial、
実行producer全bytesと旧lowparameter report bindingを保持する。
reusable runnerのisolated `--global-onecut` modeは旧十四FunctionDef ASTをすべて保持した。
gate PASS、120秒/1536MiBのfinite contract、実測2.426秒/63,328KiB。
CASは上の有限恒等式を検証したもので、独立幾何reviewは別に必要である。

残る三fixed source（ss common K零、nil T-square K零、nil common K零）の
三/四cut義務はこの一cut証明で代用しない。
