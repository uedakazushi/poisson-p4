# 二HF全族は、完全有限graphで得られるsmooth24親へ入る

2026-11-08、PRIMARY `E_reduced_nonnormal`。producer算術と全formal-support証明。
別作者literal graph再構成は未実施。二親構成の同一性／既MPS名／新成分加算は未決である。

## 入力と、計算前に証明した有限性

`HF_EQUIVARIANT_FINITE_WHOLE_SLICE_PROOF.md` の全証明を使う。
同じ二つの固定anchorだけでfullD101/Euler24/M48、unique normal character−2を得ていた。
Aのnormalized torus入力character範囲[−8,6]、Bは[−6,4]。
旧101minorのcoordinateUとhomogeneous normalNで真の102sliceを作り、
selected101式のformalIFT graphの一意性とC*equivarianceから、
**全completionのgraphそのもの**がA degree四／B degree三のliteral polynomialになる。
新しいjet次数、profile、point、minorを計算結果に応じて追加していない。

第二rank二planeへの直接既知bridgeも先に検査した。
AのSingDは唯一plane、Bはplaneとtwisted-cubic coneであり第二planeがない。
両planeのnormal位数は二なので、既wholeconcurrent m三定理をそのまま流用しなかった。

## 旧minorと真のnormal

Mの既存48minor inverseで、既存D-kernel49のcharacter−2部分からNを一つ選んだ。
全51actual group columnsを殺すlinear functional、N上のexact非零値を保存した。
AのNは元kernel index13でfunctional値75/16、Bはindex14で値2。
全DN零、N character−2、旧101/48determinantの再検証を行った。
新rank/minor profileを作らず、同じcomplete kernelのbasis quotientを読んだだけである。

## 全literal graphの結果

同じ101implicit block inverseで、Aではj2,3,4の三段、Bではj2,3の二段を一回ずつ解いた。
各補正はU内でcharacter−2jだけに支えを持つ。
full graphの150coefficient×degree配列、全補正、normal、全weightsを保存した。

| source | 全graph次数 | 元全350式の全検査次数 | 全residual |
|---|---:|---:|---|
| HF k−1/2,p2/5,ε0 | 四 | 零から八 | 恒等零 |
| HF k−2,p−1/2,ε0 | 三 | 零から六 | 恒等零 |

selected101式も全次数で恒等零、残249式を一つも落とさない。
元350式がPi係数にquadraticなので、上表以上のt次数は存在しない。
weight boundと一意IFTにより、未確認higher formal termも存在しない。
従ってこれは有限jetの推測ではなく、**全formal sliceがA1であることの完全証明**である。

実groupの48directionsと102sliceのproduct mapは、原点で微分が可逆なのでétale。
Poisson schemeの局所completionは48smooth group directions×C[[t]]となる。
rawはsmooth unique49、Euler24とprojective scale一を除いて、
各anchorは**smooth unique24dimのfull projective component**上にある。
この段階で親の名前又は二親が相異なるとは主張しない。

## 実曲線と全元familyへの戻し

Aはt=144e、Bはt=3eで全graph係数が整数になる。
この最小weighted denominator-clear scaleと全e一の十括弧をreportに保存した。
全eについてliteral actual Poisson曲線であり、e零は元の同じanchorである。
normalizedC*はt→s^(−2)tを与えるので、e非零の全点は同じactual projective PGL orbitとなる。
étale productのdimension下界とこのtorus作用から、そのorbitは24dimで親componentにdense。
generic e一の新しいtangent matrixを追加する必要はない。

旧whole image上界23とirreducible PGL×A10 sourceの閉image論証は、
各元HF全A10の像閉包がanchorの23orbitclosureと一致することを既に示している。
親componentはclosed/PGL invariantでanchorを含むため、その全orbitclosure、従って全元A10 familyを含む。
全source parameter零／境界を含めて戻し、familyの一般fiberだけに結論を限定しない。

八HFの残六源はfactor／元square-rowで空だった。
従ってR2L/m二constantnormaldoubleの残HF源は、この二actual24親の全点へ被覆された。
既知componentへの最終認識は、e非零のdense24actual点の実PGL同定又は実SNC退化などの別証明を要する。
新親を無断でMPS40又は新componentとして数えない。

## 実行identity

数学的成功はNEW V2 reportだけである。
V1はmetadataのelse23 spacingによるSyntaxErrorでcompile時に停止し、数学算術を開始していない。
failedcaller/job/receiptは不変に保存し、成功又は追加jet実行と呼ばない。
V2は同じpoints/minors/degreebounds、metadata spaceとidentityだけを直し、事前compileとgate後に一回実行した。
2.973秒/RSS52456KiB、全PASS。
shared6022/helper sourceは編集せず、exactinlinecallerをarchive/job/reportへ保持した。
