# tangent conic・m2 double root：一回の全source判定と残る一源

2026-11-07。主義務G。producerの有限QQ/QQ(t)算術。
全tangent-leading模型への幾何的還元は別pureproof依存であり、独立review前のconstant-root正常形を
本算術だけで認証しない。全E帰属及びnormal-double全Gは未証明。

## 1. 対象と全方程式

`PLANE_M2_DOUBLE_ROOT_CONSTANT_NORMAL_FORMS.md` 第2、3節の三fixed模型、及び
`PLANE_M2_MOVING_DOUBLE_ROOT_NORMAL_FORM.md` 第5節の全t nonzero族だけを一回で処理した。
triangle/concurrentを追加profileとして混ぜない。
三fixedのraw fieldsはN=2y partial_x+z partial_y、f=z(xz-y^2)に対し

| 模型 | X | Y | Z | k20 |
|---|---|---|---|---|
| h=f,Z0 | -2N | -2N | 0 | 1 |
| h=f,ZN | -2N | -2N | N | 1 |
| h=z^3 | -8N | -5N | z partial_x | 0 |

各々のactualG2は3f u^2、3f u^2、-9z^3u^2と全150から直接確認した。
k20 nonzeroという一律のunit仮定は使わない。
moving族のactualG2は3t z(yu+zv)^2で、t=0だけが別m三以上の補集合。
全P0で元350 Jacobiが零。

全weight一30、weight二9の既忠実Euler商源を使い、全350行のJ1,J2を再構成した。
full39 source variablesの全quadratic projective Jacobiと全actualDもuniversal polarizationsとして保存し、
linear-boundに不要なcompatibilityを捨てたsourceとは呼ばない。
新job `jobs/plane_m2_double_tangent.json` は先にgate PASS、実算術は2.943秒でPASS。
reportは `reports/CLASSIFICATION_PLANE_M2_DOUBLE_TANGENT_20261107.json`。

## 2. actual filtered normalizerからprojective像上界へ

Levi groupは(GL3 x GL2)/C*、dimension12。
P0のprojective stabilizer dimensionをsとする。
全25 Euler columnsのrankは24であり、Levi全13 columns、さらにP0 scaleを加えた全150行matrixGでは

\[
 s=37-\operatorname{rank}G.
\]

各matrixのcomplete kernel、元全行での消滅、rankとnonzero minorをreportに保存した。
scalar GL5中心を二重に数えない。

Sを保つfiltered groupは更に全六base-from-normal shearsを持つ。
この六次元groupはP0をactual tensorとして固定するとは限らないが、associated-graded P0を固定し、
全P1/P2 source内で作用する。忠実30/9 Euler商へ戻す操作は線形gaugeなので同じfiberを保つ。
従ってG=PGL5とfullhigher fiberから全Poisson空間への写像のfiberにはdimension s+6のgroup graphが入る。
任意のmodel parameter pを戻すと、全fullsourceのprojective像は

\[
 24-(s+6)+p+\dim\ker J_1+\dim\ker J_2\tag{1}
\]

以下である。
P1は全ker J1に入り、固定P1上のP2はweight二Jacobiのinhomogeneous affine fiberで、
空又はker J2のtorsorの部分集合だから、このdimension上界は非線形Jacobiを解かずに有効。
全higher compatibilityは像を小さくする方向にしか働かない。
成分性、smoothness又は全scheme同型は主張していない。

## 3. 三fixed模型の結果

| 模型 | ker J1 | ker J2 | s+6 | 式(1) | 結論 |
|---|---:|---:|---:|---:|---|
| h=f,Z0 | 8 | 0 | 10 | 22 | 全像22以下 |
| h=f,ZN | 7 | 0 | 9 | 22 | 更に全actualDがu因子、指定既約源は空 |
| h=z^3 | 9 | 3 | 9 | 27 | 未閉 |

第二模型では全complete K1のbad u-row v^2 coefficientが零。
全P0 u-rowはu倍、全P2にはnormal columnが無いので、全fullsolution u-rowがu倍になる。
従って実Dはuで割れ、既約性と矛盾する。
これはselected slice又はgeneric pointのfactorではなく、全higher solutionへの必要因子である。

第三模型の27は粗上界であり、27次元成分の存在でもGの閉鎖でもない。
全12higher必要parameter及び全39 quadratic Jacobiが保存されている。

## 4. moving族は全complex tで像23以下

generic QQ(t)ではker J1=6、ker J2=0、s+6=9、parameter dimension一を含め像22以下。
選んだ三rank minorのnumerators/denominators及びcomplete K1/K2の全denominatorsをQQ[t]で因数分解し、
t=0を除く全irreducible壁を元matrixの専門化で処理した。
全壁は次の七rational値だけであり、非有理complex rootを無断で落としていない。

| t | ker J1 | ker J2 | 像上界 |
|---|---:|---:|---:|
| 6 | 6 | 0 | 21 |
| -3/2 | 7 | 0 | 22 |
| -5/2 | 6 | 0 | 21 |
| -3 | 8 | 0 | 23 |
| -9/2 | 6 | 2 | 23 |
| -6 | 7 | 1 | 23 |
| -9 | 7 | 0 | 22 |

全壁でs+6=9、壁自身のparameter dimensionは零。
従ってこの**whole moving-root sourceの像は23以下**。
seedのT18所属からhigherへ飛躍する推論は使わない。

## 5. 保存と次の数学

共通runnerはisolated `--m2-double-tangent` modeだけを追記し、全旧11 FunctionDef ASTが不変。
直前のexact source b45036...及び同runnerを参照する全五current契約bytesは
`inputs/plane_m2_double_tangent_frozen_inputs_20261107.json` に保存した。
旧report/source/contractのimmutable bytesは変えていない。
current旧jobのcompatible source binding更新は統合担当だけが行う。

これで指定四源中三源はG非優越又は空へ閉じ、一源h=z^3だけが未閉。
契約通り次はこの唯一12parameter必要sourceのOre又は別の構造数学へ移る。
同じ全profile、minor列、GB又は高次jetを反復しない。
全tangent-leading削減と本算術・normalizer橋の独立reviewが揃うまでglobal claimへ採用しない。
