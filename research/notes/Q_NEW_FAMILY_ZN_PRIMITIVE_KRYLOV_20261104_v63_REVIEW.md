# ZN のdegree10原始性：有限free模型の独立検証

2026-11-04。対象は `Q_NEW_FAMILY_ZN_PRIMITIVE_KRYLOV_20261104_v63`。
同じ保存済みoctic-base点index0と同じprime32003だけを用いる。新点・新primeは使わない。

元の16 βと受理済みquadratic C/D/moment tensorsから値を独立に計算し、pencilは原著の
24 permutationとは別に、traceless4×4のNewton式
`det(C+tD)=tr((C+tD)²)²/8−tr((C+tD)⁴)/4` で全五係数を求める。
F4とp0のprime-unitを検査し、R,Z0,Z1,Z2、Gの全891項、保存済みC1の全21,012項を
独立にspecializeする。G leading coefficientは単元であり、monic化したringはbasis
`1,r,…,r9` のfree10である。

accepted ZN式から、まず小さいs多項式

```text
N=Z0*r*(r−1)²*J*NE1+25*Z0²*Fn*NE2,
D=20*r*(r−1)²*Fn*J²
```

を全係数で展開する。両者をdegree4で直接homogenizeし、sの分子は `−C1_constant`、
分母は `C1_linear` とする。原著のfactorized JH/FH/NH計算とは別の順でU,Vを再構築し、
全100 matrix成分を比較する。special algebraでC1_linearやVを逆転しない。
行列式は別のfraction-free Bareiss condensationで確認する。

## 正のcharacteristic-zero証明

baseの環を `Z_(32003)[R,Z0,Z1,Z2]` の保存residue点で局所化し、G leadingを逆転したものと
する。monic Gはfree10代数を与える。U,Vもその代数の元で、列
`V9,U V8,…,U9` の係数行列はこのbase環に定義される。
residue determinantの非零は、この**係数函数**がcharacteristic zeroでも非零であることを
証明する。空性や次元のmod-pからの持上げではない。

generic algebraは既受理の既約Gによるdegree10のfield。元のchartとC1回収からVはそのfieldで
非零である。従ってgenericにV9を除けば列は `1,ZN,…,ZN9` の共通非零倍となり、独立性から
ZNの最小次数は正確に10。**specialized Vのunit性は不要**である。
従って `Kshape=QQ(R,Z0,Z1,Z2)(ZN)`。

full NEW target体との一致と被覆次数は別のapplication packetで、全16源の斉次性・
degree14非零不変式を使って証明する。ここでは新しい不変式展開やNEW定義式を作らない。
