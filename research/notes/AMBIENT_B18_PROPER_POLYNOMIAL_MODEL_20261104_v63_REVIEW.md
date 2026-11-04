# B18 proper-source14/17模型の独立監査

対象は PROPER_POLYNOMIAL_MODEL_PREFLIGHT と PROPER_POLYNOMIAL_MODEL の凍結済み二 packet。
作者の builder や source helper を import せず、physical ey-wedge vector field の式から
CA/CB 十九基底を別構成する。divergence がそれぞれ -f_h, -2f_h であることから Euler
補正を直接適用し、全150係数、両 source center、accepted P0 を確認する。

GL4 は (h,x,z,w) を動かし y を固定する。全16方向の physical action の rank15、
parabolic11方向の rank10、actual source-gauge15 minor -2/3、同じ physical zero cut15
の minor32を検査する。GL4 中心は -P0 を動かすため、affine scalar は既に含まれる。
projectivization を追加せず、十五変数の group section だけを使う。

smooth proper zero sources の各 dimension24 と、P0 上の同じ normalization source への
同一視は既存の独立受理済み前提を保持する。source gauge15 の unit Jacobian により
全48 source germ は group section15 と gauge slice33 の形式積となる。incidence ideal
は GL4 不変なので slice から引き戻される。同じ group15 に対する physical cut15 の
unit Jacobian で群変数を形式的に解くと、actual zero cut と source gauge は全 scheme
として形式同型になる。incidence の滑らかさや reducedness はここで仮定しない。

CB19 の定数左逆と完全 annihilator17を検証し、全360 transformed36 coefficients と
全150 image/readback identities、全170係数多項式と14変数17式を再生する。原点の
Jacobian rank5・minor8/3も元の literal equations から求める。

Npair は I=(v1,a+v0,b,c,d,e)。全17式への exact substitution は零であり、五選択式と
実 flag d が六 conormal classes を張ることを非零6 minorで検証する。Nakayama により
completed local ring で I=(選択五式,d)。従って五 unit IFT の後の残る十二式は、この
実 d を因子に持つ。作者が別途作業中の rational elimination、graph membership、CM、
socle、D0 depth はこの監査からは結論しない。対象は normalized tensor E4 であり D0
との同一視を行わない。

新しい source jet、GB、local elimination、数値点は実行しない。実行結果・時間・入力
hash は report.json と independent-proof.json に保存する。
