# 非被約R・m四generic211：全実Dを二固定＋一parameterへ

2026-11-07、主義務G。producer手証明と有限linear Casimir課題。
Rのsupportが二直線以下である全generic211枝を対象とする。
新しいfull-source rank結果は後続finite jobの実測を待つ。

## 1. 全leading spaceとcollision支持

rank二plane S⊂SingD、actual Dは被約既約quintic、Sに沿う位数は四とする。
generic leading binary quarticのroot partitionを2+1+1、coefficient image rankを二以上とする。
`PLANE_M4_REPEATED_LEADING_GEOMETRY.md`§§2–4の全P1 Bertini及びsymmetric-pencil論証を使う。
これらの節はR被約性を使わず、parallel G4/fのroot partitionがU=S\Rで変わらないことだけを使う。
normal fixed double rootをuとしG4=u²Q、Q=A u²+2Buv+C v²とすると、
CとT=AC−B²の両collision divisorのsupportはSupp Rに含まれる。
unused componentsは許す。

rank(A,B,C)三ならTはirreducible conicであり、Rが2L+M又は3Lという非被約cubicに反する。
rank二のregular symmetric pencilはsemisimple又はJordanの全二型。
semisimple型でmarked double rootが両eigenvectorsと異なれば三つのused linesが必要なので同様に不可能。
残る全leading formsはactual normal/base GLの下で

`sem: G4=u⁴A+u²v²B`、

`Jordan: G4=u²v²A+u³vB`

の二つである。A,Bはindependent base linear forms、係数の全nonzero amplitudesはこの選択に含める。
どちらも二つの相異なるused linesを要するので、R=3Lの全high-rank211枝はこの時点で空である。
R=2L+Mではどちらのlineが二重かを固定せず、両orderingを残す。
image rank一及びgeneric22/4は既essential-three全八sourceへ、generic31は別の全-R theoremへ送る。

## 2. semisimple collision：全actual Dの二正常形

actual D=u⁴A+u²v²B+g5(u,v)と書く。
Aへnormal linear formを加えるactual shearでg5のu5,u4vを消す。
Bへの同じshearでu3v2,u2v3を消せる。
従って残る全sourceは

`u⁴A+u²v²B+αuv⁴+κv⁵`。

既約性はκ非零を強制する。κ零なら全Dがuで割れるためである。
normal diagonal scales、A/Bの補償scalesと全体尺度によりκ一とし、
残るu/vの自由な比でα零又は非零を0/1へ戻せる。
全actual Dは

`Fε=u⁴A+u²v²B+εuv⁴+v⁵, ε∈{0,1}`

に尽くされる。全逆はminus shearとnonzero diagonal inverse。
vへuを足すchangeでleading pencilを保つとは仮定しない。

F0はactual linear field Y=−4A∂A−2B∂B+u∂uを持ち、Y(F0)零。
したがってそのprojective PGL4 orbitは次元十四以下。
F1についてlinear stabilizerの存在は必要なく、orbit次元十五以下を用いる。

## 3. Jordan collision：一つのwholeparameterは残る

actual D=u²v²A+u³vB+g5とする。
A/Bへのnormal linear shearsがu4v,u3v2,u2v3の全係数を消す。
残る項はγu5+αuv4+κv5である。
γ零なら全Dがvで、κ零なら全Dがuで割れるので、両方非零が必要。
normal diagonalと全体尺度でγ=κ=1へ戻すが、その後のratio自由度はfinite fifth rootsだけとなる。
従ってαを一律0/1にする過剰normalizationは使えない。
全actual Dを

`Fτ=u²v²A+u³vB+u⁵+τuv⁴+v⁵, τ∈C`

という**whole complex line**として保持する。
元のcoefficient invariantはα5/(γ κ4)であり、τのfinite ambiguityは像次元を増やさない。
全τでactual Ga field Z=u∂A−v∂BがFτを保存するので、fixedτのPGL4 orbitは十四以下。
τ零も含め、specialτの増大stabilizer又は増大Poisson fiberを捨てない。

全三displayed familyはAに関してprimitive linearなのでC上絶対既約。
四gradientのconstant independenceからessential-variable数は四であり、第五座標wは全vertex方向である。
この確認はcandidate actual Dの幾何であって、Poisson存在又は全点成分帰属の主張ではない。

## 4. 必要最弱のfull linear Casimir sourceとG target

canonical divergencezero strict quadratic liftのPfaffian–Casimir identityから、
固定D class [F]を持つ全projective Poisson sourceは

`P(LF), LF={quadratic Π in five variables: div Π=0, Π#dF=0}`

に入る。全150quadratic coefficientsを保持し、vertex w dependenceを落とさない。
div0はcanonical representativeの正当な条件であり、任意Euler代表をglobal stabilizerとしたshortcutではない。
identity自体は`ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md`§1に証明され、essential-fourにも同じ証明が使える。
full nonlinear JacobiとactualPf=Fを課すことはこの必要linear sourceを狭めるだけである。

essential-four baseの選択Gr(4,5)は四次元。
sem F0の全D familyは十四＋四で十八以下、sem F1は十九以下。
Jordan wholeτは十四＋四＋一で十九以下、fixedτでは十八以下。
従って必要なLF dimensionsと全Poisson像上界は次の通り。

| 全actual D | 十分なlinear source上界 | 全projective像 |
|---|---:|---:|
| sem F0 | dim LF≤7 | ≤18+6=24 |
| sem F1 | dim LF≤6 | ≤19+5=24 |
| Jordan genericτ | dim LF≤6 | ≤19+5=24 |
| Jordan finite specialτ | dim LF≤7 | ≤18+6=24 |

一回のfinite jobで全五成分・degree六のCasimir係数1050とdivの全linear係数25を作る。
固定二族はQQ、JordanはQQ[τ]の全1075×150 matrixとし、complete generic kernel/minorを保存する。
generic minorの全irreducible factorsを尽くし、必要なら各QQ[τ]/pの全exact characteristic-zero fieldでrankを認証する。
finite wallの全conjugatesは一つのfield certificateで覆い、浮動小数又は一部のrootだけを使わない。
得られた必要源が大きい場合は全sourceとfalse criterionを保存し、actualPf/rank四geometryへrouteを変える。
full30/9 connection-profile反復、GB、長いjet又はpoint farmはこのtargetには不要。

## 5. 射程と残余

本noteの全normal formsはR被約性を使わないので、その実Dを持つ全Poisson sourceに適用出来る。
ただしnonredR高rank211の全幾何からG結論を出すには、上記finite source dimension判定の成功が必要。
m二/三、m四squarefree、generic31及びimage-rank一は別packetへ分けたままである。
zero cubic Rの全支持又は全分類の完了は主張しない。
