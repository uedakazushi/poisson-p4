# 22個の小minor多項式のQQ4局所化イデアル、一回のGB

対象は前段で実sourceとの等号を検査した22 generator（2963項）と、49項の
U=adFBCJによるzU−1だけ。ringはQQ[a,e,d,z]、全次数1、GRevLex。既存M2 1.24.11の
default gbを一回だけ使い、ChangeMatrixを保存する。過去に問題があった
AlgorithmLinearAlgebraは指定せず、環境の再セットアップもしない。

unitなら同じcached GBから1//Gで元23 generatorのBézout liftを取り、M*L=1を
検査する。保存したQQ4係数を別のnative FLINT積和でも検証する。properなら
全basisとcached source-change witness M*T=Bを保存し、同じnative checkerで
元generatorへの所属を確認する。independent GB recipeは固定basisの全S-pairと
元23 generatorの剰余検査であり、GBの再計算や点による推測を要しない。

60秒・2GB・32MiBでinput/log/phaseを必ず残し、時間境界や未完了liftをunit証明と
しない。B rank3の判定であって、rawH7やfull239の自動証明ではない。properの場合も
primary decompositionや新しい点の選択はしない。
