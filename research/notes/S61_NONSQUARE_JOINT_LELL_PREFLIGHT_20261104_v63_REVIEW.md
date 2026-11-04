# joint L/ell input と98 unit pivots の独立監査

2026-11-04。対象は `S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63`。
四 quad を単項式 ideal と同一視する案は受理しない。full source の Z/W 補正を含む
actual six quadrics と、その monic leading terms による正確な Schur recipe を監査する。
新しい rank、Hilbert function、GB、core点、標数、genuine/quartic 全係数展開は使わない。

## 1. source と座標を独立に復元

元六 quad は既存の actual sparse source から作る。physical変数 v を

    (v0,v1,v2,v3,v4,v5,v6,v7)
      =(x0,x1,x2,W,x3,H,Z+x1−3x0,x4)

へ反復 sparse polynomial multiplication で変換する。ell は実六 quad の恒等式
ell W=−(de+8)q1+8q2+dq3 の W係数として回収し、保存 ell7 と全係数を比較する。
L は Z=W=0制限の q1 derivative により、4C*Lj=−405*(∂q1/∂x1)j を検証する。
これは full q1=L²/8100 を仮定する操作ではない。

new Y=(L,ell,x0,x3,r,Z,W)、r=7x2/2+3x4 の forward matrix F をこの元 source から構成する。
保存49 entriesと照合し、保存 inverse numerator N について FN=NF=delta I7 を全98係数で検証する。
また、この matrixの determinantは四つの coordinate rowを除いた固定3×3式

    L1*(3ell2−7ell4/2)−ell1*(3L2−7L4/2)
      =delta=3200*a*C*d²/81

そのものであることを確認する。従って N/delta は FULLU 上の actual inverseであり、
未知 minor の inversion はない。

## 2. full six quadrics と非零 overlap

新六式は p1=8100q1、p2=8100(e q1−q3)/(400a)、p3=ell Z、p4=ell W、
p5=q0、p6=q5。後二つの ell Z/W source identity も元六 quadから係数ごとに確認する。
q0=p5、q1=p1/8100、q3=(e p1−400a p2)/8100、q5=p6 の後に
q2=(p4+(de+8)q1−dq3)/8、q4=((10a+3e²−3e)q1−3(e−1)q3−ap3)/(10a)
と戻せるので、row change は既知 a と有理数だけの可逆変換である。

レビューは author の forward substitutionを使わず、保存 p numeratorの全174係数に
**newY=F*oldX と逆向き代入**し、元 q0,…,q5 と上記定義へ戻ることを検証する。
保存した各 row denominatorは非零有理数と a,d,Cだけの積であることを厳密に確認する。
結果は1633 core terms・72 nonzero entries。

最初の四行の leading terms は descending lex で L²,Lx0,ell Z,ell W、いずれも monic。
しかし p1,p2 は計82項の非零 Z/W補正を持つ。さらに四つの単項式に属さない degree2 monomialの
非零係数も存在するので、この fixed frame で actual idealをその四 monomial idealに置き換えることはできない。
S=x0p1−Lp2 は保存 numerator/denominatorと全91係数で一致し、72項・10非零係数の非零式。
一方 Wp3−Zp4 は全91係数で零。二つの cubic overlapを同じように零とは扱わない。

## 3. unit triangular square と全 overlap の保持

weighted monomialsは重み (1,1,1,1,1,1,1,2) の指数再帰から独立に生成する。
四 prime quad と degree(D−2) cofactorの全積を作り、その distinct leading monomialごとに
先の source rowを選ぶ。保存された26行/98行、target pivot列、complement列と全て比較する。
各選択行の全 supportについて、leading coefficientだけが対応する pivotにあり、
他の項はそれより lex小さいことを確認する。従って選択 squareは対角1の三角行列で、
inverseは有限和 Σ(I−U)^k。これに未知 determinant factor は含まれない。

degree3は target91→65、four-quad products28→overlap2。
二つの overlapは source row7=p2*L と row26=p4*Z。
後者だけは選択行 p3*W と literalに等しく、その Schur residualが零なので削除できる。
q0/q5倍行14、genuine16を保つと正確な31→65の source-span recipeになる。
これは元の product41 basisと同じ source spanに、既知の一つの零関係を復元した表示からの操作である。

degree4は four products116から98 unit source/target pivotを除き、**全18 overlapを保持**する。
さらに q0/q5 products58、元 genuine16の old-light倍行112、全33 true quarticを保持する。
従って actual full319→target239 は、同じ cokernelの **source221→target141** に移る。
ここで genuine112とquartic33の source labelを変えず、quad174だけを新 generator/cofactor basisへ変換してよい。
targetの light coordinate substitution は全行へ共通に作用する。
cofactor basis changeを genuineの old-light multiplier labelsへ誤って適用しない。

65/141は unit squareの complement数であって、actual quotientの Hilbert functionではない。
all18 overlapを保持したため、単項式割り算で nonzero S等を失うこともない。

## scope

本 REVIEW は actual入力・既知 unit frame・26/98 pivot recipe と全 source/target countを受理する。
genuine16や33quarticの新しい core展開、221×141係数行列、rank141、full239、全点cubic57、
free183 follow-onはこの監査で構成・証明していない。後続183表示は別の有限 circuit監査に分ける。

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_joint_Lell_preflight_review_20261104_v63.py --output work/marked-square-v63/joint-Lell-review-REPLAY
```
