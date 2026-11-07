# Ga9零c三軌道：固定nil22 unfoldingの不適合と幾何的制約

2026-11-07。`E_reduced_nonnormal`。
全五源の有限化は `NIL211_SHORT_TRANSLATION_REVIEW.md` で独立受理済み。
ここはそのうちc=0三源の特定変形経路の診断であり、成分帰属・非所属を判定するものではない。
各全十一源の閉包は保存されたcanonical A2=w²点の23次元PGL軌道閉包に等しいので、
その三点を固定して実変形を探した。

## 1. 指定したcofactorの一回のunfolding

元の `N0=y∂x+∂z` に対し

    Nε=y∂x+(1+εw)∂z

とする。ε≠0では元nil22-image座標を

    (X,Y,Z,W)=(x,y,z/ε,w+1/ε)

と取れば `Y∂X+W∂Z` の実affineGL移送である。
全27係数nil22源を保持し、各係数のLaurent次数−4,…,2の線形消去を一回行うと、
325×189係数系の階数126に対し三targetともaugmented rank126であった。
従ってregularな**線形源**極限はある。しかし自由係数を零とした最短候補は全てJacobiを満たさない。
この段階ではPoisson族も既知成分への包含も得られていない。

## 2. 全cotangent形のactual pencilは存在しない

次に一つの完全な実pencil条件 `Pε=P0+εR` を解いた。
base=(y,w)、fiber=(x,z)とし、Rの四anchorを全二次式

    Rxy=U1(y,w), Rxw=U2(y,w), Rzy=V1(y,w), Rzw=V2(y,w)

とする。全24係数を使い、`Ryw=0`。
divR=w∂zを保つ全fiber-fiber式は

    Rxz=R0(y,w)+(w+divV)x−(divU)z

であり、R0の六係数も全て許す。
cotangent二次元baseのJacobiではR0は関与せず、全cross/self Jacobiでその六微分が零であることも確認した。
RがPoissonなら `[divR,R]=0` なので、追加の必要線形式
`−w divU−U2=0` を使える。
全cross Jacobiとこの必要式を解き、その完全解へ全self Jacobiを戻すと次になる。

|canonical源|全線形系rank / augmented rank|残るanchor自由度|本質的障害|
|---|---:|---:|---|
|quadratic_w_power1|23 / 23|1|self Jacobi係数に非零定数63/8|
|quadratic_w_power3|22 / 23|なし|cross Jacobiの時点で不整合|
|linear_w|22 / 22|2|self Jacobi係数に非零定数2|

従ってこの指定Nεとf={y,w}=0を保つ、ε一次の**全actual cotangent pencil**は三点とも存在しない。
全24行列係数、全cross Jacobi、全必要線形式、全self Jacobi、power3の線形左functionalと
非零pairingは `CLASSIFICATION_GA9_NIL22_UNFOLDING_DIAGNOSTIC_20261107.json` に保存した。
focused checkerは `code/classification_ga9_nil22_unfolding.py`、実行0.664秒。
有限契約 `jobs/ga9_nil22_unfolding.json` はgate済み。

これは同じ固定cofactorに対する高次実変形、別のcofactor変形、既知成分への帰属を排除しない。
線形源候補のJacobi失敗とactual pencil失敗を受け、この経路は停止した。
より長いLaurent列や高次Jacobi jetを自動的に追加しない。

## 3. 通常の固定normalconeでも三軌道は得られない

受理済みの全ambient軌道階数48から、各canonical点のPGL軌道次元は
`48−Euler24−scalar1=23`。従ってPGL安定化群の次元は1である。
一方、非零冪零linear場 `Nhom=y∂x+h∂z` が二ベクトルを保存するので、
その安定化Lie代数はちょうどC·Nhomである。
この一変数Lie代数には非零半単純元がなく、安定化群に非自明torusは含まれない。

任意の非自明PGL一径数部分群λ(t)によるprojective軌道極限は、そのλ(C*)に固定される。
従って、この三canonical軌道は**どの固定tensorの非自明PGL1-PS極限にもならない**。
固定linear centerに沿う通常normalcone、又は対角weight極限だけを探す経路は初めから不適合である。
これはMPS成分の軌道閉包への非所属を意味しない。
moving-basis arcのような非1PSの実変形や、本質的plane-motion磁気class生成なら別の問題である。

この三つの零 c 源の E 帰属は依然未解決である。非零 c 二型の解決は末尾を参照。
この診断自体による新成分は追加しない。

## 4. 固定 unipotent 一径数軌道極限も全有限次数で不適合

1PS の代わりに、同じ安定化冪零場 `N=y∂x+h∂z` による

    ε^k exp(N/ε)_*Q

を考えた。これは高次 jet の切断ではなく、表現上の L_N が冪零なので初めから有限な実共役族である。
全150係数上で Euler 射影 C を保持し、`rank C=126`、`CL_N=L_NC` と
projective 表現での冪零指数5を正確に確認した。
非零 endpoint が P0 なら、最大冪 k に対して `L_N^k Q` が P0 に比例する必要がある。

完全な有限次数の像判定は次の通りである。

|Ga9 source|P0 が入り得る正次数 k|追加の必要条件|
|---|---|---|
|quadratic_w_power1|1のみ|`L_NQ=P0`, 全350 `D_P0 Q=0` は rank111/112 で不整合|
|quadratic_w_power3|1のみ|同じく rank111/112 で不整合|
|linear_w|1又は2|下記により両者とも不適合|

linear_w の k=1 では全線形解の実自由度は18である。
その**完全18変数 Q** を全350 Jacobi に戻すと、パラメータに無関係な非零定数48が現れ、
Poisson Q は存在しない。k=2 では最上の cross 条件

    L_N²Q=P0,  D_P0(L_NQ)=0

も必要である。全500×150の係数系の rank は63、augmented rank は64であり不整合となる。
全ての負結果の非零左 functional 又は literal Jacobi 定数係数を保存した。

証拠は `CLASSIFICATION_GA9_UNIPOTENT_LIMITS_20261107.json`（3.140秒）、
`CLASSIFICATION_GA9_UNIPOTENT_LINEAR_W_20261107.json`（1.085秒）、
`CLASSIFICATION_GA9_UNIPOTENT_QUADRATIC_LIMIT_20261107.json`（1.974秒）である。
対応する三つの `classification_ga9_unipotent_*.py` と有限 job は全て gate 済み。
全 nilpotence 範囲を尽くしたので、この固定 unipotent 共役経路はここで停止する。

これは三源の MPS 帰属を排除しない。Ga9 三源の安定化群の Lie 代数は C·N なので、
別の純 unipotent 一径数部分群を同じ endpoint で選ぶことも同じテストに帰着する。
従って単一の乗法型又は加法型一径数部分群の固定 tensor 軌道極限では足りず、
一般の moving-basis arc 又は別の全源構造が必要である。
なお Ga9 の非零定数 c と二根 c の全源帰属は、それぞれ
`GA9_CONSTANT_COMPONENT.md` と `GA9_DISTINCT_QUINTIC_COMPONENT.md` で別に解決済みである。
