# 八つの平面自由五次源：新七候補の全周囲成分性・独立監査

2026-11-07、reviewer `root`。全seven ambient算術 **PASS**。
曲線リストの尽くしと全源の認識は別作者の独立監査を待つ。
このnote又は計算reportだけで新しい成分数を登録していない。

`classification_first_essential_quintic_review.py --eight-ambient` は、
既公表の全九FunctionDef ASTを不変にした独立modeである。
producer、legacy helper、SymPyをimportせず、既存の別作者Fraction疎演算だけを使った。
有限180秒/1536MiBの事前gate後、5.542秒・30224KiBで全証明書を再構成した。
旧104番目のC45は凍結済み全ambient/二次支持証明を再利用し、不要な再計算はしなかった。

## 1. 全七点と全150入力

全候補ごとにSaito二vector fieldから磁気零・GL2単位frameのpointを別構成し、
integer scale後の全150係数をproducer anchorと照合した。
全350projective Jacobiを再計算し、さらに各150basis方向との和から全350×150微分Dを復元。
Euler二十五・actualGL5二十五・scalar一の全51列Mも別に計算した。
DM=0、M rank48の非零minor、Dの全非零minorを検査した。

| 曲線型 | rank D | projective tangent | 全支持の最弱証明 |
|---|---:|---:|---|
| A12 | 102 | 23 | tangent下界一致、smooth23 |
| E8+A4 | 101 | 24 | 一normal、二次非零 |
| E6+A6 | 101 | 24 | 一normal、二次非零 |
| E6+A4+A2 | 101 | 24 | 一normal、二次非零 |
| A8+A4 | 102 | 23 | tangent下界一致、smooth23 |
| 3A4 | 102 | 23 | tangent下界一致、smooth23 |
| A6+3A2 | 102 | 23 | tangent下界一致、smooth23 |

各Dのkernelは表示48列とnormal零又は一列で全て張られる。
全input/output character、Dのintertwining、actual48columnのweightと独立性も確認した。

## 2. genuine48cutと必要な二次証明だけ

全七点でactual48basisの48coordinate minorが非零、normalはその48coordinate上零。
implicitな101/102標準座標、全48actual directions、normalからなる全150直和のdeterminantを
別に計算し、真のlinear cutであることを確認した。

三つの一normal点では、選択101式の101×101逆行列について両積を検査した。
normalの全350二次Jacobi、選択IFTの一意二次補正、全350残差を再構成し、
全150列Dを消す余核functionalと非零pairingを確認した。値は順に

`−47546093/13602816`, `−913230352584/7808230609235`, `−3431736/3246473`。

一変数formal cutは非零n²初項を持つので、被約支持は原点。
raw local上界48とactualorbit/gauge/scalar下界48が一致する。
よって各23orbit閉包は全周囲の23成分となる。
四smooth点もraw tangent48と同じ下界でsmooth23・局所唯一成分となる。
三次、正確な冪零指数又は完成環全体の計算は不要で、実行していない。

## 3. 相違のscopeと残る独立入力

各actual quinticのgradient span rank三を直接確認した。
immutable prior104はcommit `6ab28a327b64543ba0a7bf8cd302bea9def517e6`、
catalogue blob `f7ec76f1b2bb1850f2e3314e36f673f218558c1c` の全bytesで照合した。
prior103の相違監査は、既受理の全MPS actualSNC witness、旧二既約五次のessential5 minor、
他のrank/dimension/closed factor条件が、今回の全essential3・既約・23成分に一様に適用する。
その独立入力のbytesを固定し、旧29MPS計算を理由なく再実行しない。

旧104番目の最大四重cuspに対し、新七候補は全て最大multiplicity三以下。
八つ相互の28pairは異なるcusp配置を持ち、essential三変数の内在的base quotientにより
PGL5同型ならPGL3曲線同型となる。二つの稠密23orbitは同じ23閉包を持てない。

この最後の曲線型認識と、全必要曲線リスト及びfullGL2+六磁気源の独立完全性は別監査の入力である。
そのPASS前は新七つを採録しない。全G/E・成分数上界・文献新規性も未証明のまま保持する。
