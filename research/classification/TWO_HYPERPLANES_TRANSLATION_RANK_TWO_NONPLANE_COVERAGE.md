# B=0・rank二平行移動：非平面二全源のE被覆

2026-11-08、PRIMARY `E_reduced_nonnormal`。
全source有限化、全normalizer、22軌道閉包橋、二つのactual半単純曲線を
別作者のFraction算術と幾何読解が一括してPASSした。
新しいambient成分は加えない。

## 1. 正確な対象と有限化

実退化因子Dが相異なる二単純超平面と一つの既約三次因子からなる型1+1+3を考える。
二Hをh,xとした忠実全模型で、x-normal fieldは

```
V=(B y+c x+d)∂y.
```

本定理はB=0、c,dが独立で、Dが非平面である全支持を扱う。
可逆fiber変更によりV=x∂u+h∂vとなる。
初めからWの低次項を零とせず、元の全45二次係数から線形Jacobiを解いた全21係数sourceを使う。
入口は `TWO_HYPERPLANES_CUBIC_SUPPORT.md` §§1,5.1と
`CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_20261107.json` のj111_rank2である。

A=n11,B'=n18と書くと、三次の非平面部分は(Ah−B'x)(hu−xv)。
(A,B')≠(0,0)の全支持は可逆normalizerと全Jacobiのunit消去により次の二源へ尽くされる。

| source | 完全parameter域 | 三次のlinear正常形 |
|---|---|---|
| L=h | A8のn17非零開 | h²u−hxv+x³、本質四変数 |
| L=h−x | 全A7 | (h−x)(hu−xv)+hxw、本質五変数 |

全係数の必要十分性と全逆は
`TWO_HYPERPLANES_CUBIC_TRANSLATION_RANK_TWO_NONPLANE_SOURCES.md`。
A=0,B'≠0はH/x交換・Euler gauge・tensor−1で第一源へ移る。
両非零はw尺度とx/u尺度で第二源へ移る。
n17=0は第一源の三次にh因子を生むため既約三次の対象から除外される。
A=B'=0は平面sourceであり、本定理では未処理の補集合として保持する。

## 2. 全sourceはそれぞれ同じ22軌道の閉包

`TWO_HYPERPLANES_TRANSLATION_NONPLANE_ORBIT_CLOSURES.md` の証明を使う。
rank2・square-zero cofactorのprojective class[N]は11次元で、
P(im N^dual)=P1内の相異なる二H covectorsはさらに2次元を持つ。
固定exactNの非平面Jacobi fibreは各9次元なので、各全sourceのPGL飽和像は22次元以下。
同じ構造的kappa=1点で実Euler/GL/scalar像rank47、Euler鉛直rank24から実projective orbitは22次元。
全source微分もrank47で、各source及びPGL5は既約である。
従って各source飽和閉包Sは厳密にその22軌道閉包となる。

この橋は全parameterを任意の一つのorbitへ送る同型を仮定しない。
稠密orbitの閉包が全Sに等しいという次元と既約性の論証で、全境界も含む。
実ambient D96・七余法線はこのpositive membershipに不要であり、
七法線の新しい二次map又は高次jetは実行していない。

## 3. 二つの真の既知親と同じcanonical点

同じ二canonical点は、別々の全三parameterの実半単純sourceからの曲線で得られる。

| target | semisimple cofactor | 全parameterと実曲線 |
|---|---|---|
| L=h−x、Perazzo | (2,1,−1)、j1,ell0 | TWO_HYPERPLANES_TRANSLATION_PERAZZO_SS_CONFLUENCE.md |
| L=h、本質四 | (3,1,−1)、j1 | TWO_HYPERPLANES_TRANSLATION_ESSENTIAL_FOUR_SS_CONFLUENCE.md |

各親は全parameterでstrict Poisson、同じ既約A3 source内に実SNC rank四点を含む。
その全像を含む一つのambient成分はSNC点にも会うから、
全SNC成分尽くし定理によりMPS40の一つである。
したがって親の全像が既知MPS有限閉和へ入る。

表示された真のLaurent parameterをそのまま親へ代入し、inverse shearとtensor epsilon倍、
固定PGLを施すと全負epsilon冪が消える。
全kappaについて元の全150端点、strict10・元350 Jacobi、全Pfaffian、全逆を確認した。
L=hでは既に指定されていたcanonical PGLの逆も用いて、元の同じ点を全150係数で復元した。
source parameterの積を独立化せず、kappaも一度も反転しない。

epsilon≠0は既知親の真の点、epsilon=0は正則な同canonical点である。
閉性からkappa=1が既知閉和に入り、§2の全orbit閉包橋によりS全体が同閉和へ入る。
kappa=0は同じliteral曲線で直接含まれ、その他のparameter境界はSの閉性で含まれる。

## 4. 独立確認と結論

`TWO_HYPERPLANES_TRANSLATION_NONPLANE_SEPARATE_REVIEW.md` と
`CLASSIFICATION_TWO_HYPERPLANES_TRANSLATION_NONPLANE_SEPARATE_REVIEW_20261108.json` は、
producerを実行又はimportせず、全21source/全Jacobi係数、全消去とnormalizer逆、
実orbit/source rank47、13+9上界、両全三parameter親・全epsilon/kappa・150/350/Pf/SNC/全逆を確認した。
一つの有限Fraction batch、0.233秒、RSS23932 KiBでPASS。

従って**B=0かつrank(c,d)=2の非平面全1+1+3支持はMPS40有限閉和へ被覆される**。
二源閉包の次元は22であり、既知MPS成分の最小次元23より小さい。
よってこれらを新しい22次元ambient成分として数えることはない。

平面A=B'=0、他のnilpotent B、既存p=0源、Ga9 linear-w、stable二源、short r2/3及び全Eは別に残す。
旧globalC39/C40陰性と七法線の停止記録は保持するが、このpositive全源被覆の仮定ではない。
