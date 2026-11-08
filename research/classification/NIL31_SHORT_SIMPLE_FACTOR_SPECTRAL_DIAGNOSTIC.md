# nil31短鎖 r=2/3：simple-H companion の whole endpoint 障害

2026-11-08、PRIMARY E。root提案の新 simple-factor 入口を独立に全係数確認し、
一回の有限 exact source matching を行った。
**この fixed block-spectral frame は target に達しない。全r2/3 EはOPEN。**
因子を選び直す入口は正しく、旧 multiple-H の negative proofs とは別の幾何である。
今回の fixed frame を任意の coupled spectralflag の normalform とする定理は無い。

## 1. 成功した正確な入口確認

`NIL31_SHORT_SIMPLE_FACTOR_SPECTRAL_ENTRANCE.md` の P−E∧L、scale−3/7、
H=旧Z,x=旧X+15t旧W/7,y=旧Y,z=旧W,w=旧h を t=d1で独立構成した。
全150係数で新 target と一致、Z-row零、元 projective Jacobi零も確認した。
これらは raw report の `simple_factor_gauge_transfer_all150_exact` と全target vectorに保存した。
新 infinity の quartic は w²(wy−Hz) の倍数で、不可約ではない。
既存 nil22long irreducible-quartic 七源への無条件包含は使わない。

## 2. one bounded run と保存した全源

Nε=y∂x+(H+4ε²x)∂y+w∂z+ε²z∂w の generic eigenvalues は(-2,-1,1,2)ε。
全 spectralframe は

    X=H+4ε²x−2εy,Y=w−εz,Z=w+εz,W=H+4ε²x+2εy。

det−32ε⁴、explicit inverse は入口ノートと公開inputに固定した。
全SS12九係数、全七 genuine product relations、全63 Laurent coefficients[-4,2]を保持した。
endpoint/no-negative-polesの必要条件は150×63の exact QQ系となり不整合。
rankは係数51、拡大52である。
この段階で nonlinear relation solver は実行せず、actualcurve 350/Pf検査にも進んでいない。
gate後の一回の実測は0.154013秒。

保存matrixには次の unit row functionals がある。

| row | 旧新target側の coefficient | 全63column | RHS |
|---:|---|---|---:|
|12|Pxy[wy]|零|27/49|
|21|Pxy[Hw]|零|−3/7|
|22|Pxy[Hz]|零|−27/49|
|35|Pxz[w²]|零|−3/7|
|36|Pxz[zw]|零|−9/49|
|66|Pxw[w²]|零|−9/49|

row12一本で t!=0 の whole canonical target を排除するのに十分。
全matrix/RHS/RREF/constant-frame全10basisを raw reportに保存した。

## 3. 任意meromorphic source coefficientsにも通る最弱手証明

旧新座標の ∂x/∂y成分を持つのは spectral ∂X/∂Wだけである。
従って新 source を pullbackした Pxy は、全 nine parametersについて正確に

    Pxy=[n1H²+n2(w²−ε²z²)
         +(2−n7−n8){(H+4ε²x)²−4ε²y²}]/(16ε²)。

これは全 SS12 source の PXW と explicit inverseから従う。
Pxy には wy,Hw,Hz が全て欠ける。
targetの wy coefficientは27t/49で、t!=0なので矛盾。
parent n_i を任意 Laurent/rational/formal meromorphic functionsへ増やしても欠落は同じ。
旧 endpoint Laurentwindow又は選んだ nonlinear sliceによる偽の非存在結論ではない。
ただしこの block-split frame の外へは一切外挿しない。

## 4. 未解決の幾何と保存範囲

long(X,W)/short(Y,Z) spectral blocksの間をcoupleするmoving flags、
singular affine centers、変動spectralratio、別knowncomponent sourceは残る。
単にこのblock frameのparameter次数・negativepole予算・higherjetを増やしても失われたwyは出ない。
したがって同frameの再試行はしない。
正の全source被覆には、実在する別parentのwholecoefficients又は新flag構造が必要である。
新成分、generic SS12閉包からの非所属、全r2/3 E完成は未主張。
