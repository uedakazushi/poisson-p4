# rank二quadric cone：全Hamiltonian pencilのgraded fiber上界

2026-11-07。主義務G、producer plane_collision_review。以下は有限全source計算の契約であり、
実行結果の全共鳴を確認する前にG≤24を主張しない。

## 1. 最高symbolから正当に得るOre形

全sourceは `QUADRIC_CONE_RANK_TWO_HAMILTONIAN_PENCIL.md` とその全350証明書に従う。
全四parameter（radial trace sとF_vの三係数）上の最高空間次数Jacobiの全rowspaceは
`CLASSIFICATION_QUADRIC_CONE_RANK_TWO_HIGHEST_SYMBOL_20261107.json` に保存した。
その式からsemisimple CではF_v=κu1、nilpotent CではF_v=κu2が必要十分、
C0ではF_vは自由となる。全三型でconstant pencil Dは零である。
従ってこれは入口で誤ってD0を課す操作ではなく、元の全Jacobiから得た帰結である。
以後全条件は正しく

`Jac(B)=0, Lie_T B=0, Lie_Y B=Y∧T`

となる。F_v条件以外の全低次parameterと全Rの境界を保持する。

## 2. 三角形分解

semisimple/nilpotent Cの全sourceはbase(κ,s)の上で

```
B=B2+B1+B0,    T=T1+T0,    Y=Y3+Y2+Y1+Y0,
Y3=-u q.
```

B2,T1はbaseだけに依存する。全既存sourceの自由座標を用い、
weight1の十四座標は(B1,Y2,T0)、weight2の十二座標は(B0,Y1)、
weight3の三座標はY0である。parameterの一つが複数項に現れる場合も同じweightを持ち、
部分源へ制限していない。

空間次数dのJacB係数はparameter weight3−d、Lie_TBは2−d、
Lie_YB−Y∧Tは4−dでhomogeneousである。
各weight jの方程式の、weight jの未知係数に対する部分だけを取ると、
係数はbaseのみのmatrix Mjになる。
他の項はweight<jの既決parameterだけに依存する。
従って固定base上の**全**Poisson fiberのdimensionは

`nullity(M1)+nullity(M2)+nullity(M3)`

以下である。各段階はaffine-linearな必要条件であり、右辺compatibilityや残るweight4以上の条件は
fiberを空又は小さくするだけである。全sourceがこれらの必要条件を満たすことを
保存済み全original350と全weight identityから確認する。

## 3. C0の全F_v補集合と実normalizer

C0かつF_v0では、全actualPfaffianのh0係数は恒等零になる。
従ってactualDはhで割れ、既約quinticの仮定で排除する。
この主張を全sourceのactual126係数で確認し、旧rankzero C0除外を流用しない。

C0,F_v非零はSO(q)×非零scalarの全二軌道、F_v=u1又はu2で覆う。
これはq-dual norm非零と零の全補集合である。
実座標u'=rGu,h'=th,v'=dvとtensor scalar λでA(h,v)=1を保つ条件は
λtd/r²=1であり、F_vはλGによって動く。任意のλをtで実現できるので、
F_vの大きさを一にする操作は全source上可逆である。
全Fのu-only五係数と全hLも可逆に動き、zero divisorの種類を制限しない。
両軌道のaffine dimensionsは三と二、固定後のbaseはs一係数である。
同じ十四・十二・三の全未知を保持して各一変数operatorを使う。

AutS×tensor scalarは十三次元、Aのsem/nil/C0 strataは七/六/四次元。
Asem/nilを固定した全sourceのdimension dに対するPGL像上界は18+d/17+d。
C0からさらにF_vをsem/nilに固定するとstabilizerは六/七次元であり、
同じ18+d/17+dとなる。dにはκ,s又はsの全base dimensionを含める。
これは全actual automorphismのfree frame fiberによる上界である。

## 4. 全共鳴を残す有限証明書

各Mjのgeneric非零minorと完全kernelをQQ(κ,s)又はQQ(s)上で保存する。
二変数minorの全因子が線型なら全lineをQQ(tau)へ引き戻し、そのgeneric minorの全複素根を
有理数又はexact number fieldで再評価する。lineの交点も全てこの有限coverに含まれる。
非線型因子又は時間上限に達した場合はそのまま保存して停止し、point/profileの追加を行わない。
base dimensionをgeneric二・line一・point零へ下げた上でnullity和を足す。
一変数C0もgeneric一・全根零を全て保持する。
全strataのPGL像上界が二十四以下と確認出来た場合だけ、whole rank二coneのG非優越を結論する。
