# unit98 Schur と actual free183 表示の独立監査

2026-11-04。対象は凍結済み `S61_NONSQUARE_JOINT_LELL_FREE183_20261104_v63`。
joint L/ell の actual source と monic98 squareは直前の independent PREFLIGHT REVIEWで受理済み。
以下を typed scalar/linear expression checkerで独立に構成し、保存された全 operator rootと照合する。
author script、rank solver、GB、追加点、183×141係数全展開は使わない。

## 1. source174の両方向変換と target NF

新 targetの light 座標は (L,ell,x0,x3,r,Z,W)、Hは不変。
FとN/deltaから physical oldv8↔newY8 の linear mapsを構成する。
Z=3v0−v1+v6 を fullに含める。degree2の全29 cofactor monomialの像は、二つの linear forms
の積または単独Hとして構成する。old→new、新→oldの全 monomial-map rootsを照合する。

qprime=Cqold の六 source change と明記した逆行列を、両積がI6になる scalar identitiesで確認する。
cofactorの変換は contragredientであり、

    cprime_j = Σ_i (C^{-1})_{ij} cold_i(oldv(newY)),
    cold_j   = Σ_i C_{ij} cprime_i(newY(oldv))

となる。old174→prime174 と逆向きの全 linear rootsをこの規則から構成する。
source basis changeは **quadratic174だけ**。genuine112の multiplierは元の old v_j、
全33 quarticの source labelも元のままとし、その polynomial target coordinatesだけを変える。

六 prime quadから全174 product coefficient rowsを構成し、98 selected source rows、98 target列、
残り76 source rowsと141 target列を照合する。全18 overlapを残す。
U98、RD98×141、NC76×98 の全 rootを照合し、Uが対角1の上三角であることも確認する。
rowの target projectionは μU=row_C を前進再帰で解く

    μ_j = row_C,j − Σ_(i<j) U_ij μ_i,
    NF(row) = row_D − μ RD

である。保存された239 independent formal target入力→141 outputの全 linear rootsが一致する。
これは未知 minorを割らない exact Schurなので、full319の cokernelを221→141へ移す。

## 2. restorationとPnew/Knew

prime quadratic sourceの残り76係数 x を入れ、selected98 sourceの係数を betaとする。
pivot targetが零になる条件は U^t beta+NC^t x=0。従って

    beta=−U^{-t}NC^t x.

レビューではこの符号と転置を含む再帰、prime174への挿入、old174への逆変換、
その後の pinned P22適用を全 rootで照合する。これが Pnew:76→22。
旧 sparse K22:22→174を prime sourceへ変換して76へ制限した Knewも全 rootを照合する。

Pnew Knew=I22は新しい巨大な22×22展開をせず、次の exact operator proofによる。
旧Kは pure-quadratic sourceの literal kernelなので、変換したKも全239 target係数が零。
unitUにより、これのselected98係数は残り76から上のrestorationで一意に復元される。
従って old sourceへ戻すと旧Kそのものになり、受理済み P22 K22=I22が適用できる。

**この復元の同一性を使うのは pure-quadratic K22 上だけ**。
一般の full kernelはnonquadratic rowsの pivot-target寄与を持ちうるので、その二次部分だけを
restorationして元の二次部分へ戻るとは仮定しない。この仮定はPnewの構成にも補空間の証明にも不要。

## 3. 既知38方向とsource183

quad174だけを変換したため、旧D15の quartic11/cubic4 pivot座標と全ての nonquadratic
source coefficientsはそのままである。unit98 Schurはsource kernel上ではselected98係数を
除いて制限するだけで、unitUによりkernelを同型に移す。従って同じ15 unit source deletionsが使える。
残るsourceは quad76とsurplus130。

独立受理済みfree281 packetのold15処理後の新vectorについて、quad174だけをprimeへ変換し76へ制限する。
u22,h108は以前と同じliteral nonquadratic roots。
Pnew/Knewによる

    k=Pnew*x76−u22,
    z76=x76−Knew*k,
    kernel184=(z76,h108)

の全 scalar rootsを照合する。これはPnew Knew=I22により、76+22+108 sourceの可逆分解から得られる。
自由184表示は [A76+Bnew Pnew,Cnew]。sourceが固定した元184行の部分行列だという主張ではない。

separatorはliteral Cnewの座標76,77に残る。同じfree281 scalar prefixと同じ零identity rootを使い、
その和は既知unit theta=−1024000aC³J/243。scalar prefixの一致を確認し、全追加atomの分母も
既知 a,d,F,B,C,Jだけからなることを確認する。
座標和の定数変換後、thetaで正規化したkernelの成分が1なので、一つのsource generatorを削除できる。
全正規化kernel rootsと183 retained indicesを照合する。

これで **FULLU上の actual full319 cokernelは自由183→自由141の cokernelと等しい**。
全既知方向は38、残り42のまま。rank141、full239、全kernel exhaustionは未証明。
全2705 linear nodesは純粋に線形で、入力はold174/prime174/rest76/kernel22/target239の
五つの disjoint rangeに型付けし、各operatorが予定されたrangeだけを使うことも検証する。
全大規模scalar rootはfactored expressionのままで、nativeは1012 atomの既知unit分母と
前と同じ小さなseparator identityだけを展開する。

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_joint_Lell_free183_review_20261104_v63.py --output work/marked-square-v63/joint-free183-review-REPLAY
```
