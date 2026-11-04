# 独立後検証：114 source 五次係数と全619 pure cubic値

本監査は提供された係数を代入して検証するもので、新しい未知係数の求解
ではない。作者scriptをimport/実行せず、別の Fraction series reader、
symbolic finite-ring arithmetic、実 source column再構成を用いる。

## 114係数の後検証と quotient二次値

CA/CBそれぞれ元の57未知の保存u0..4と、提供されたu5を合わせ、元の57
scalar方程式をQQ[r]/r六乗で全部評価する。各familyで全57式の全6係数が
実target（最初55が零、最後二つが既証slope*r）と一致することを確認する。
保存lower4と逆行列を持つ固定J0により、この五次係数は一意である。
独立なu5計算を重複して行う必要はない。

元の signed numerator/denominator rootと actual common r曲線への輸送は
別の入力独立監査917f249…へ固定する。global h0=0、両g4=0、CB curvature2=0
を保つので、ordinary係数の関係bbar=g*hからh2=bbar5/g3を再計算できる。
r三乗を取り除いた後の非零有理数g3だけで除し、零を反転しない。
全bbar=g*hの係数をdegree5まで再照合する。

## 全619 source functionalの独立再構成

既証の完全Schur/Theta criterion724db717…の actual Theta、Ystar2_rr、U0、
Ystar1、m0、piから、quadratic source adjointをforward Eulerで再構成する。
作者の31-entry配列は最後にだけ照合する。元620 constant source結合と
実a1、Xp1、および619 quadratic-kernel basisから、rawP55三次、Xp55二次、
二次source correctionと既知補正を全部再構成する。全619列のsource labelと
係数が凍結配列へ一致することを検査し、単に作者の619数値を読み上げない。

旧612列のf63 multiplierと、proper gA planeのbasisを区別し、literal
source nodeIDを保持する。nativeのsource55は全次数で零である。二次native
代表は零-origin結合にした後だけ使い、個別nativeX2の値と誤認しない。

## 既存source3と全未知感度の扱い

実common105の保存二変数source3をrA=sA=rへ縮約し、元105式全てをdegree3
まで確認する。actual YB7/8へownCB57 source3（保存source4の部分）を合成し、
コピーした元57方程式全てを同じ次数で確認する。新しいsource solveはない。

rawP3を評価する際、既存の全regular Hadamard/cofactor germsの定数・一次
を使い、それ以外の二次と三次を独立な形式変数として保つ。値を零とした
というsource主張はしない。ordinary係数の和、積、単位逆を通じて感度を
運ぶ。未知同士の積はdegree4以上なので、この計算では厳密に寄与しない。

二次source correction側も、B11上のtyped sum/productをforwardに読み、
同じnodeIDと同じr²方向を同じ形式変数として扱う。受理済み108 composite
source値を用いる部分はその独立source輸送へ明示的に依存する。cached値に
ないatomの二次係数を勝手に補わない。

619個の最終式で、ordinary h2/h3およびg3の全感度が消えることを確認する。
残るのは八つのprimitive g2のr²感度だけで、それも別source packetの行列と
比較する。known部分を既存cache-economyの619成分へ照合し、その615個の
missing-annihilator結合も独立に零と確認する。最後に、今回後検証した八つ
の実g2値を代入して全619値を判定する。

## 完全性と結論の範囲

全619実値が零なら、既証の完全initial-space criterionによりpure r³候補は
生きる。actual HFはdegree3まで(5,3,2,1)、全長は11–17となる。これは
一部のspecialized kernelや一部sourceだけの零判定ではない。

真の全長、さらに高い層、DVR torsion exponent、全nil三乗の零性は決めない。
既証のdepth19は保持する。元の五次producerの結果は独立監査待ちと記録
されていたが、本後検証はそのfrozen bytesを変えず別reportで受理する。

## 再現

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_T_gdiff_second_univariate5_review_20261104_v63.py
```

60秒・2GB・32MiB・150000参照のguard。新しいsource jet、未知求解、rank探索、
348 scalar行列は用いない。
