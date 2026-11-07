# Twisted cubic cone・rank二支持の全G非優越：別作者による独立検証

2026-11-07、reviewer plane_collision_review。**全PASS**。
producer/root の `CUBIC_CONE_RANK_TWO_BASE_REPRESENTATION.md` と四つの凍結reportsを独立に監査した。
新reviewer `classification_cubic_cone_rank_two_review.py` はproducerをimportしない。
幾何polynomialsは別実装のstdlib Fraction sparse arithmetic、三whole λ moduleの係数だけQQ(λ)を用いる。
GB、jet、有限体又はλのsample列は使っていない。

対象はtwisted cubic cone

`S=(x0x2−x1²,x0x3−x1x2,x1x3−x2²)⊂P4`、vertex x4方向、

を非正規支持に持つrank四Poissonで、実Dは既約被約quintic、π|Sは一般rank二である枝。
全vertex、Rの全可約・非被約を含む。
結論は全PGL5飽和射影像のdimension≤24というG非優越だけである。
Eの全点帰属、成分数及び全G分類の完成は主張しない。

## 全150入口と完全27源

全150 quadratic homogeneous coefficientsについて、各qiに対する
`xi*(Π#dq)_j−xj*(Π#dq)_i` をconeのdominant parametrization
`(s³,s²t,st²,t³,c)`へpullbackし、全690行を再構成した。
intrinsic bracketは

`H_ij=x0{qi,qj}−2qj{qi,x0}+2qi{qj,x0}`

であり、raw bracketへ置き換えていない。その全五partialsのpullbackが全285行。
元975×150の**全entryとrow labels**がreportに一致した。
cone idealはこのparametrizationのprime kernelなので、projective tangencyの必要条件の両方向が正しい。
smooth dense locusのintrinsic conormal可換性はHのfirstnormal消滅と同値であり、この全条件を使っている。
vertexに追加の局所条件が必要なら源はさらに小さくなるだけで、上界に影響しない。
vertexを除いた小sliceのimageだけを全源と呼んでいない。

全51 kernel columnsについて975行との積をQQで全て零と確認し、columnsのrank51を別Fraction RREFで確認した。
指定99×99 minorを別消去で計算するとdet=**−192**。
したがってrankは99、kernelは正確に51である。
全25 Euler wedgesを再構成してrank24と全rowとの積零を確認。
canonical divergence-zero 27columnsの全150係数、全divergence零と
`rank[ Euler24, canonical27 ]=51`を確認し、全projective源が27次元であることを両方向で得た。

## 全350 Jacobi、全126 actualD、base九＋vertex十八

27源の全parameter quadratic monomialsをpolarizeし、元の全350 projective Jacobiと全126 actual Pfaffianを再構成した。
base indicesは正確に `[0,1,2,3,4,5,6,7,16]`。
各base tensorはcに依存しない四変数Pと`c div(P)∧∂c`、残十八は
c-independent、divergence-freeな`V∧∂c`であることを全係数で確認した。
全vertex-vertex Jacobiは零。
base-baseとbase-vertex termsは元350の異なるspatial rowsに現れる。
従って全Jacobiの零条件は、base九変数の十二quadricsと全350×18のhomogeneous vertex moduleの零条件に正確に分かれる。
全Dもbase-base＋base-vertexの全126行として保存・比較した。

旧packetの「c-dependent Dがbase Jacobi rowspaceに入る」は**false**であり、
その全非零remaindersを独立に再構成してfalseのまま保持した。
この旧diagnosticを今回の上界の仮定に使っていない。

## Actual SL2と全base支持の被約尽くし

四座標のSym3 actual SL2、c固定のh/e/fをそれぞれ全150係数へ直接作用させ、
九源上の全三action matricesを確認した。
producerの明示9×9座標変換とinverseの積はidentity、変換後は三つのSym2標準作用である。
十二元quadricsの全変換を再構成し、binary quadratics U,Vの積UVの全五coefficientsが
元十二式の明示constant linear combinationsであることを確認した。

C[s,t]はdomainなので、複素被約支持でU=0又はV=0。
両whole branchesの全十二式を、次の単純式と**両方向rowspace同値**で比較した。

* V=0：U∧W=0。
* U=0：V∧W=0、およびdisc(V)=0。

非零UではW=λUとなり、Uは相異二根又は平方の二型。
非零Vは平方でW=λV。
U=V=0ではWの相異二根・平方・零を全て残す。
binary quadraticの二相異根又はdouble rootをactual SL2で移し、nonzero tensor scalarで係数を一にする。
逆は同じSL2行列の逆とnonzero scalarの逆であり、vertex十八も可逆に移る。
λの全complex lineを残しているので、この正規化は選択した点だけの検査ではない。
Wも零ならbase P=0で全Π=V∧∂c、rank≤2となりrank四枝から除かれる。

この全補集合から、rank四源は正確に三whole λ familiesと二fixed familiesに含まれる。
固定λでのbase strata dimensionではなく、全orbitとtensor scaleを含めて
A_ss=4、A_nil=3、B_nil=3、W_ss=3、W_nil=2である。
有限λ wallsではこのdimensionが一減る。

## 全五modulesと十九walls

各canonical familyの元350×18 matrixの**全entry**を、独立に得たbase-vertex coefficient rowsから構成した。
全126 particular actualD coefficientsと全126×18 vertexD matrixも再構成した。
全kernel vectorsの350行との積零、free-coordinate block identity、指定maximal minorの非零を検証した。
kernelによるrank上界とminorによるrank下界が一致するので、全moduleのrankとkernel dimensionは完全である。
generic whole λではQQ(λ)の厳密なfield operationsのみで検証し、λをsampleしていない。
minorと全kernel denominatorの全irreducible factorsを独立にfactorし、reportのwall集合との一致を確認した。
全十九factorsはlinearなので、その正確なQQ根で**元の全350行**を戻し、再び全kernel/free block/minorを確認した。
全generic/rational wallのminorは零でない。

| case | generic rank/kernel | exact walls | wall kernel最大 | 全射影像上界の最大 |
|---|---|---:|---:|---:|
| A_ss | 14/4 | 5 | 6 | 24 |
| A_nil | 14/4 | 7 | 6 | 23 |
| B_nil | 14/4 | 7 | 5 | 22 |
| W_ss | 14/4 | 0 | — | 22 |
| W_nil | 14/4 | 0 | — | 21 |

minorのzeroを単にrank dropと推測せず、すべてのwallを元matrixで再検証した。
逆にgeneric minorの非零開と全denominator開の外には、上のfactor以外のλは無い。
λの無限遠はW-only二型へ既に含まれ、base零も別処理済みである。

## 全像のdimensionと正確なreview範囲

embedded coneのPGL stabilizerは八次元：twisted cubicのPGL2三、vertex相対尺度一、
cへの四base linear formsのshears四。
全sourceはこのactual groupで保たれ、PGL5×source上の対角作用は第一因子でfreeである。
個々の点のstabilizerが大きくても、八次元fiberを引く論理は失われない。
各base stratum dimension dとkernel kから

`dim whole projective image ≤24−8+d+k−1=15+d+k`

を得る。finite wallはdを一減らす。
上表から全image≤24、これでtwisted cubic coneの全rank二支持はG非優越となる。
actualDのproper factorを使う追加の推測は必要無い。

実行17.956秒、RSS62,796KiB。
gate/hashは保存整合性、全行列・kernel・minor・functionalの独立再構成と上記whole geometryが数学的reviewである。
producer原本・旧false reports・主原稿・global catalogue・Gitは変更していない。
