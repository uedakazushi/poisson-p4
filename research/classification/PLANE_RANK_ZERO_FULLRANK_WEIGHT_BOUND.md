# rank零plane・m=2・O(1)係数rank三：全像dimension二十三以下

2026-11-07、主義務G。producer exact arithmetic及び幾何接続。
独立review前は新global claimへ混ぜない。

**命題。** `PLANE_RANK_ZERO_MULTIPLICITY_TWO_FIRST_JET.md` の全m_D=2、
common co-Higgs saturation O(1)、coefficient map M:C³→sl2がrank三の分岐は、
全projective Poisson像のdimension≤23を満たす。
全rank零plane、m_D≥3又は他のco-Higgs型を含む命題ではない。

## 1. 全normalization sourceと全higher-weight係数

前掲noteの全Lax還元と実可逆base/normal変換により、
τ≠0のsemisimple/nilpotent二whole-symbol familyに到達する。
base(a,b,c)、normal(z,w)、M=[[b,c],[a,−b]]、E_n=z∂z+w∂wとし

`Π0=T∧E_n+U∧((bz+cw)∂z+(az−bw)∂w)+K∂z∧∂w`。

| 型 | U | T | K |
|---|---|---|---|
| ss | ∂b | 2τ(a∂a−c∂c) | (1−2τ)zw |
| nil | ∂a | τ(c∂b−2b∂a) | (τ−1/2)z² |

全30係数Π1及び全9係数Π2を前掲noteのcanonical Euler quotientで保持する。
各D_iは全350 projective Jacobi coefficientのlinear operatorで、
`D1=E∧[Π0,Π1]`、`D2=E∧[Π0,Π2]`。
非線形なweight二以降の方程式を捨てても必要条件として

`dim(higher full Poisson fiber)≤dimkerD1+dimkerD2`           (1)

が成り立つ。
fixed Π1に対するΠ2方程式はD2のinhomogeneous systemだからである。
さらに高いweight Jacobiはこのsupersourceを狭めるだけである。

## 2. 正確な一parameter certificate

reusable runnerは `code/classification_plane_rank_zero_fullrank_weights.py`。
既受理runner `classification_plane_smooth_cubic_weights.py` の純粋なSchouten/Euler helperだけをimportし、
既存code/reportのbytesを変更していない。
job `jobs/plane_rank_zero_fullrank_weights.json` は全input hashと120秒/1536MiBの契約を持ち、
classification.py gate PASS後に実行した。

`reports/CLASSIFICATION_PLANE_RANK_ZERO_FULLRANK_WEIGHTS_20261107.json` は
全350×30 D1、全350×9 D2、実shear30×6、完全kernel、nonzero maximal minorを保存する。
各symbolic kernelへの全行代入はQQ(τ)上で零。
minorはQQ[τ]上でexactに求め、全irreducible factorを尽くした。
全factorはlinearであり、全rootをQQ上で再計算した。
generic minorの零点を無視したgeneric-rank主張ではない。
実行は4.871秒、最大RSS55,816KiB。

最初の同じ算術は両family・全rootまで成功後、factor multiplicityのSympy IntegerをJSONへ書く時に停止した。
runtime reportは元producer bytesを保持する。
唯一のserialization int変換を宣言・gateした後、同じ有限算術を再実行して成功reportを保存した。
算術の失敗又は異なるsourceへの点列probeではない。

## 3. generic rankと全例外根

両familyでgeneric rank(D1)=24、rank(D2)=9、actual shear rank=6。
従ってgeneric kernel dimensionsは6,0である。
以下の表は選んだrank certificate minorの**全**非零例外根である。
minor固有の余分な根も削除せず、全350行でrankを確認する。

| 型 | τ | dimkerD1 | dimkerD2 | higher upper d | 全projective image upper |
|---|---|---:|---:|---:|---:|
| ss | 1,3/10,1/4,1/6,−1/6,−3/2,3/2,−1/2 | 6 | 0 | 6 | 22 |
| ss | 1/2 | 6 | 1 | 7 | 23 |
| nil | −1/2 | 6 | 0 | 6 | 21 |
| nil | −3/2 | 7 | 0 | 7 | 22 |
| nil | 1 | 6 | 2 | 8 | 23 |
| nil | 3/2 | 8 | 0 | 8 | 23 |
| nil | 1/2 | 7 | 1 | 8 | 23 |

τ=0もcertificateでは保持し、両kernel8,0となるが、m_D=2 stratumの外であることを明示した。
これをempty又はm≥3の解決と扱わない。

generic selected-minor product（非零constantを除く）は

`ss: τ⁵(τ−1)(2τ−3)⁴(2τ−1)⁶(2τ+1)⁵(2τ+3)³(4τ−1)²(6τ−1)²(6τ+1)²(10τ−3)`、

`nil: τ⁵(τ−1)²(2τ−3)¹¹(2τ−1)¹²(2τ+1)(2τ+3)²`。

よって未検査のalgebraic exceptional rootは無い。

## 4. 全base parameterの復元

fixed planeでのfirst-symbol source dimensionsはss十二、nil十一。
前掲noteのp二・M九・τ一、nil discriminant一式という全parameter countを使う。
τを一つのrootへ固定した層ではそれぞれ十一、十となる。
plane choiceはdimension六、tensor projectivizationは一を引く。
この直接parameter sourceに(1)を足すので、generic image upperは

`ss: 12+6+6−1=23`, `nil: 11+6+6−1=22`。

固定rootではそれぞれ `16+d`、`15+d` となり、全表で二十三以下。
sourceに実現しない非線形Π1/Π2を含めたupper boundで既に十分なので、
nonlinear support分類又はall350 GBを追加しない。
finite rootsのunionとgeneric openで全symbol familyを尽くすため、命題が従う。

## 5. generic枝では実finite shearによりD因子も確認出来る

reportのshearはactual coordinate change `x→x+zY+wZ`、Y,Z∈C³の微分である。
canonical Euler quotientに直した全six-column matrixがD1 kernelへ入り、generic rank六なので全kernelを尽くす。
このcoordinate changeはunipotent linear transformationで、逆は符号を戻すだけ。
normal weightを上げ、全quadratic bivectorはweight二までしか持たないのでfinite conjugationでΠ1を消せる。
消去後のweight二方程式はD2Π2=0、generic injectivityからΠ2=0となる。

全symbolのactual Dはcertificateで

`D_ss=−4*a*c*τ*(a*z²−2*b*z*w−c*w²)`、
`D_nil=c²*τ*(a*z²−2*b*z*w−c*w²)`。

従ってgeneric shearだけの枝は既約被約quinticにならない。
この追加の空性はG≤23に必須ではなく、全special rootには前節の弱いboundを用いる。

## 6. 未処理の補集合

co-Higgs saturation O(0)、O(1)だがrankM≤2、multiplicity≥3、
rank二restrictionの別枝はこの計算の外にある。
新99成分のようにm二symbolがsquareの点もnil型の全special rootに含まれるが、
個々の成分帰属又は全E点分類をdimension upper boundから推論しない。
