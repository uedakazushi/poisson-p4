# E6 smooth-flex／A5 smooth-hyperflex：真二法線と二つの全23成分の別作者確認

2026-11-07、`/root/quartic_review`。**全算術、formal sliceの数学、全source境界及びprior124／相互相違を独立に受理する。**
新しい点、高次jet又はcompleted ringの完全表示は用いていない。
既MU6/MU5で保留していた二つの同じcanonical整数点を使う。
componenthoodは本稿の新証明であり、全marked source／曲線型は既独立証明を再利用する。
台帳・原稿・Git及び成分数は統合担当だけが更新する。

## 1. 不変入力と別実装

producerの全二法線reportは
`CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_TWO_NORMAL_OBSTRUCTION_20261107.json`
（cfedee6f…）、全依存入口は
`CLASSIFICATION_E6_A5_TWO_NORMAL_COMPONENT_DEPENDENCIES_20261107.json`（98e72a7b…）。
実行producer全文は前者に固定されており、current producer runnerはimportしない。

新checker `classification_marked_quartic_two_normal_review.py` は、既別作者の
五変数sparse Fraction演算とLaurent parserだけを使い、producer・legacy engine・CASをimportしない。
有限job `jobs/marked_quartic_two_normal_separate_review.json` のgate PASS後、
**2.430秒、22,328KiB**で両点の全算術を再構成した。
新report `CLASSIFICATION_MARKED_QUARTIC_TWO_NORMAL_SEPARATE_REVIEW_20261107.json`
に実行reviewer全文、全input hash、各unitと結果を保持した。

旧MU6/MU5全sourceから全六magnetic係数を保って同じ点を再評価し、旧integer scaleと全150係数が一致した。
両方について次を確認した。

* 元の全350×150 Jacobi微分D、旧100minorのexact一致。
* Euler全25列、GL5全25列、scalar一列。Euler rank24、全actual action rank48。
* 全六magnetic微分もこのF48に入り、全D積が零。
* 実F48の48pivot行を零にする真cut、全二法線、全100implicit方向を再構成。
  `F48+N2+U100=raw150` のunit及びD*Ni零を全係数で確認。
* 元の350式の(u²,uv,v²)三列を全展開。uv列はQ(N1+N2)−Q(N1)−Q(N2)で、cross係数を半分にしない。
* 100implicit像を消した**全350余核行**を再構成。
  選択basis formsへ全行が属す係数恒等式、及び全150D列を殺すfunctionalsからbasis formsを得る逆向き恒等式を確認。

D rankは単に大きいminorの下界ではない。
F48とN2の五十独立kernel方向、全150直和unitが上界100を与え、100minorと一致する。
全51action列と全6mag列がF48へ戻ることも48blockのexact inverseで確認した。

## 2. 群作用による真cutから全ambientへ

raw150にはEuler垂直加算、GL5、scalarのconnectedな実algebraic groupが作用する。
Euler加算はP→P+E∧L、Lはlinear。斉次二次Pには[E,P]零、[E,L]零なので
`E∧[P+E∧L,P+E∧L]=E∧[P,P]` が恒等的に成り立つ。
従って単にtangent vectorsがJacobi kernelへ入るだけでなく、全projective Jacobi idealが保存される。
GL5とscalarも全idealを保存する。

Fの48pivot行をfrとし、切片Cをその48座標の変動零で定める。
F_frはunitだからraw150=F⊕C。
48個の独立な実作用方向を持つ局所群parameterとCのaction mapは、点でinvertible ambient differentialを持つ。
formal inverse function theoremにより、全raw Poisson germはsmooth48-factorと**真cutのPoisson germ**の積になる。
群作用が全idealを保存することが、ここで形式的なtangent分解だけに留まらない理由である。

Niはcut内の全二法線。
Dの100独立標準列e_jをFへ射影して
`U_j=e_j−F(F_fr)^−1(e_j)_fr` と置く。
全U_jもcut内、DUは元のD100列と同じで、非零100blockを持つ。
従ってcutは真にN2+U100で尽くされる。
保存した150detは、Uをe_jへ戻すF-column操作と50×50補行detで独立に再確認した。
100implicit変数を解くformal graphはNi方向の一次項が零なので、全残式の二次leading termsは正確に

`Qbar=Q−D_cols*(D_rows,cols)^−1*Q_rows`

となる。これは全350式から求めたmapであり、selected subsetだけのslice assertionではない。

## 3. 全二法線の支持は原点

**E6 smooth-flex。** 全Qbarのquadratic spanはrank三。
三basis formsからu²、v²をconstant結合で得る保存identityを独立に確認した。
従ってquadratic initial idealの根基は(u,v)。

**A5 smooth-hyperflex。** 全spanは正確に

`g1=(5/2)uv`, `g2=(54/25)u²−(186624/7225)v²`

の二次元で尽くされる。全degree三Macaulay行列のdetは **−2519424/7225**。
更に

`u³=(6912/1445)v*g1+(25/54)u*g2`、
`v³=(289/8640)u*g1−(7225/186624)v*g2`

を全係数で確認した。これはv一chartだけの共通根判定ではなく、無限遠を含む全P1上の空性である。

両方で全formal残式のinitial idealがm-primaryなquadratic idealを含む。
従って真cutの被約dimensionは零、全raw local dimensionは48以下。
一方実orbitがdimension48を与える。
connected作用は各全既約成分を保つので、点を通る全componentはそのorbit閉包を含む。
局所上界とorbit下界が一致し、**実projective23軌道閉包が全ambientの既約成分**である。
Euler24及びscalar一を除いた23という次元であり、S61 sliceの次元を数えていない。
両tangent dimensionはraw50なので、scheme smoothness又は精密nil厚みの分類は主張しない。

## 4. 全六magnetic sourceと零invariant境界

全sourceのSaito kernel・allmagnetic normalization・actual inverseは、
`PLANE_HYPERPLANE_QUARTIC_MU6_SOURCE_SEPARATE_REVIEW.md` とMU5版で既独立受理済み。
今回も同じ全sourceから点と全六微分を直接再構成した。

各固定marked pairの六parameter源はirreducible affine sourceであり、
既normalizerが非零magnetic invariantのdense openを一つの23orbitへ移す。
source mapは零invariantでもregularで、全境界の像はそのorbit閉包に入る。
従って今回の全componentは元の完全六magnetic sourceとそのPGL飽和閉包に等しい。
零invariant点をgeneric pointと同一orbitだとは言っていない。

両actualDはreduced H+絶対既約quartic、essential三であることが既source式と今回のactual Euler-Pfにも一致する。
essential≤3及びlinear factor存在はclosed条件なので、その全componentにも成り立つ。

## 5. immutable prior124の全監査と二新成分相互

priorは公開commit **8b46f29ddbb9f539abc7a026e492ee25bf8ad699**、COMPONENTS blob
**8b00b8c01a3c76dde6053b918e5eaa318adc1e26** の全124行projectionへ固定した。
current台帳の将来更新と混ぜない。
全124のID一致と、次の全dispositionを独立に確認した。

| 比較理由 | 行数 |
|---|---:|
| 全component dimension相違 |39|
| generic rank二 |6|
| same23 MPSのactual SNC essential五点 |29|
| ESSENTIAL5／Ga9constantの独立再構成済みactual essential五点 |2|
| proper closed factor locus又は旧絶対既約quintic actual witness |32|
| MINIMAL_Mのdense actual orbit cofactor rank一 |1|
| A0／DPMのD5／3A2 quartic profile |2|
| 既μ6九＋μ5四のdense marked orbit profile |13|

MPS29はgeneric divisorの予想ではなく、各accepted全sourceに属するactual nonzero SNC点を使う。
D=κhxyzw、κ非零はessential五なので、新全componentのclosed essential≤3と両立しない。
旧ESSENTIAL5/Ga9constant、A0/DPM/MINIMALのprinted witnessは旧MU6独立receiptのactual計算を再利用する。
その証明書は改変・再実行しない。
旧irreducible quintic actual点は新全componentのclosed linear-factor locusの外にあり、
他の旧factor型はproper multiplication imageの閉性又はnonreduced因子の閉性で区別する。

marked十三源は既独立に各全source閉包が23orbit閉包と確認されている。
E6 smooth-bitangentだけ同じquarticなので、H上contact partition **(2,2)** と今回の **(3,1)** を使う。
その他十二はunderlying quartic intrinsic profileが異なる。
A5には既十三源で同じquartic profileは無い。
今回二つ相互はE6 unibranch (3,4) cusp、mu六と、A5二smooth枝contact三、mu五で相違する。

ここではmarked contactが境界でclosedかどうかを仮定しない。
もし同じ23componentなら、両23actual orbitsは同一irreducible成分内のdense open orbitであり、交わる必要がある。
異なるintrinsic marked profilesはこの交わりを不可能にする。
unique linear Hとunique irreducible quarticの復元はactualDの因数分解から内在的である。

## 6. 証明状態

今回受理するのは、E6 smooth-flexとA5 smooth-hyperflexの二つの全23component、
両元complete magnetic sourceの全境界被覆、prior124全て及び相互との相違である。
既MU6/MU5receiptの旧「成分性未知」は当時のscopeのまま保持し、新証明を本稿へ分離した。
他のGa9源、double normal rootのhigher source、一般G/E又は総成分数の上界は結論しない。
