# 必要五二次係数だけを保持する実 B11 source 商

実 X2 の adjoint に必要な二次 monomial は
`r²,z0*r,z1*r,z4*r,z5*r` に限る。そこで実 CA-F9 座標から

```
B=QQ[z0,...,z6,r]/(r³,all zi*zj,all zi*r²,z2,z3,z6)
```

へ写す。基底は `1,r,r²,z0,z1,z4,z5,z0r,z1r,z4r,z5r` の11個である。
これは全八変数の三次 jet ではない。

既に証明した全次数 conductor pair は、実単位変換後
`IA=(z7,r³)`, `IB=(z7+c*r²,r³)` である。全 cross product は上記 B で零と
なる。従って実 `H20=Rbar/(IA*IB)` から B への全射がある。本 builder は
H20 の20基底の全400積を独立な B の monomial 乗法と照合する。B基底の
H20 lift は `[0,8,16,1,2,5,6,9,10,13,14]` である。

元の105未知根および164 source readback の保存済み H20 値を、この環写像
で B へ運ぶ。新しい105求解はない。両 actual monic `nuA,nuB` も同じ写像で
運び、`C_B=B[tauA,tauB]/(tauA²-nuA,tauB²-nuB)` の44次元乗法を保存する。
全44三乗組85184の結合則を確認する。nuの全pair積はこの実 B で零である。
単位逆の幾何級数には保守的な最大 ideal 九乗零 guard を使える。

prospective gauge job の元の定式化は15 constraints in
`B[etaA]/etaA²=nuA`、すなわち30本の scalar B 方程式である。既存の30 scalar
gauge 函数の定数・一次を使い、五つの二次係数だけを補正するので未知数は
**150 QQ** である。補正は `m_B²` に属し、`m_B*m_B²=0` より既存の origin
J30 逆だけで足りる。新しい Jacobian 二次求解は不要である。

解いた同じ scalar coefficient を C_B へ埋め、補助 etaA を physical tauA
として使う。補助変数と physical tauA を gauge 方程式の構成時から同一視
しない。30方程式を C4 まで冗長に base extend する600 QQ定式化も可能だが、
元の scalar formulation には不要である。

将来の一回の bounded job は、残差の定数・一次が既存値と一致することを
確認し、150 gauge 二次係数、実 proper kappa/G8、native S2 を求めるものに
限る。既証の L2 cancellation により R12/T/N の二次 adapter は計算しない。
全 source ancestry は凍結 native128 と同じ actual48/105 roots へ固定する。
本 packet 自体はその gauge job をまだ実行しない。

live/saved guard は、source DAG peak3600×11、105/164根、300 trace係数、
複数の44係数 tensor/matrix scratch、固定 J30逆、600の保守的補正上界を
全部加えて150,000 QQ未満とする。出力32MiB・時間60秒・memory2GB内の
preflightであり、新しい physical source jet、旧 r/g の再除算、348行列、
T長さ・深さ判定はしない。
