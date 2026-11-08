# 自由offdiagonal connectionのsign訂正とgeneric triple全coverageの撤回

2026-11-08。明示status：**旧generic triple finite listの全coverageを撤回する。**
旧proof/report bytesは保存し、正しいminus curvatureによる全source再導出を別V2へ記録する。
全generic triple/G/Eが未解決であるという主状態は変えない。

## 正しい式と反例

normal rowsのmatrix β=[[p,0],[r,q]]はnormal polynomialにβij n_j∂niとして作用する。
このlinear-field作用はmatrix commutatorのanti-representationなので、
全ordinary-flat条件は

`dp=dq=0, dr+(p−q)∧r=0`、すなわち `dr−(q−p)∧r=0`。

旧自由r分類で使った `dr+(q−p)∧r=0` は誤り。

既検証済みmoving-double sourceでL=c/a、r=−dL+L(q−p)とすると
直接微分でdr=(q−p)∧rとなる。
具体的にf=a²b、p=(dlogf)/3、q=dlogbなら
r=−dc/a+c da/(3a²)+2c db/(3ab)。
その(a,c) curvature係数はdr=2/(3a²)、(q−p)∧r=2/(3a²)。
正しいdifferenceは零、旧plusは4/(3a²)で非零。
このP0の元full350Jacobiはmoving-double certificateで零と検証済みであり、
signを選ぶ余地はない。

## 撤回するcoverageと保持する証明

以下のgeneric triple notesの**自由r有限listが全first symbolsを覆う**という矢印を撤回する。

- PLANE_NONREDUCED_M3_TRIPLE_ROOT_FIRST_CONNECTIONS.md
- PLANE_NONREDUCED_M3_TRIPLE_FINITE_NORMALFORMS.md
- PLANE_NONREDUCED_M3_TRIPLE_CANONICAL_REDUCTION.md
- PLANE_NONREDUCED_M3_TRIPLE_PARTIAL_BOUND.md
- PLANE_NONREDUCED_M3_TRIPLE_ORE_ENTRANCE.md

誤ったlistを入力した旧90source report及び旧14Ore reportは、
**入力された各concrete firstsymbolに対する元方程式/K/Jacobi/30/9/D3/Ore算術**を保存する。
fullP0不整合/部分源の空性/部分源の上界/全Qsource恒等式は、その元入力の範囲で有効。
しかし「全generic tripleを覆う」「残余は28labelsだけ」「全未閉は14sourcesだけ」は不成立。
旧reportのarithmetic PASSを全classification PASSへ読み替えない。

影響を受けないもの：

- generic21全四primitive flags。rはparallel tensorから一意に決まり、
  r=−dL+L(q−p)という正しい式と全actual P0Jacobiで処理済み。
- constant-root obstructionと全weight-one Jacobi identity。
- PLANE_M3_CANONICAL_TRANSVERSE_OFFDIAGONAL_OBSTRUCTION.md のcanonical strict-lift全補題。
- PLANE_NONREDUCED_M3_TRIPLE_IRREGULAR_C_EXCLUSION.md の全δ非零排除。
  これは自由r-flat式を使わず、全rを残した必要式だけを使う。
- R3Lでq(C)がinjectiveな二whole potential、及びφ零/r=d(xy)の全排除。

generic tripleに関する受理済み全G閉鎖claimは存在せず、現段階の台帳全Gは未証明のままである。
このerratumは台帳を勝手に変更しない。統合担当がCURRENT/NEXT/HANDOFFの誤coverageを訂正する。

## 正しい有限入口への次手

R2Lではx=b/a,y=c/a、η=4/3−e−α、
r=(Pdx+Qdy)/x、全global P/Q八係数の正しい式は

`x(Q_x−P_y)−(η+1)Q−γxQ+δxP=0`。

log-only sourceのkernel-jump ηは−1,2であり、旧1,−2ではない。
base-normalizer resonance η一、scalar η零も別に保つ。
R3Lでは正しい式はdr−dφ∧r零。
例えばφ=xのc-dependent全r型はy dx−dyであり、旧y dx+dyではない。

正しい全r8、全normalizer/zero-fiberを別V2で再構成し、
旧actual P0が一致するものは旧算術を明示同値で再利用する。
不足sourceだけを新finite jobで検証し、旧90/14を黙って上書き又は全再計算しない。
