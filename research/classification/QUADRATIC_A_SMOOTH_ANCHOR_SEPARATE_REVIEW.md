# 同じquadratic_A λ7点のsmooth24・唯一MPS成分：独立review

2026-11-07。E_multiple_hyperplane。**PASS**。
対象は既存のλ=7、m=(1,0,0,0,0,0)、整数代表2Pという**同じanchor**である。
新しい点を選ばず、既受理の全λ実曲線の一次方向一つだけを足す。
この点でPoisson schemeはsmooth、唯一の全空間成分はprojective dimension24で、既知MPS40の一つである。
T番号は同定しない。新成分、他のλのsmooth性、別sourceの包含は主張しない。

## 1. 全sourceと曲線方向の一致

`QUADRATIC_A_SYMPLECTIC_SEPARATE_REVIEW.md` で独立再構成した全七sourceを使い、
公開元sourceのuniversal整数代表2Pと全記号のまま照合した。
λ7、m0=1、他の磁気零へ評価すると、旧weight-frame reportの全150係数と一致する。

全七微分も構成した。I=m0+m4を用いた以前の独立算術から戻すとき、
∂m0=∂I、∂m4=∂I+∂(新m4)というchain ruleを保った。
λ微分はQQ(λ) coefficientを厳密に微分してから7へ評価した。
選択した磁気方向だけに制限していない。

全λ実ε曲線をC=(q1,x,y,y−q2)の逆で元座標へ戻し、整数倍二を取った。
ε=0が旧点そのものであることと、ε微分Rの全150係数を照合した。
全350行のD R=0も独立展開した。

## 2. 全ambient tangentとactual49のexact再構成

独立Fraction sparse arithmeticで、二次raw liftの全150方向から**全350×150**のDを再構成した。
旧weight九blockの全nonzero minor及び保存inverseの左右積を照合した。
各rankは3,8,13,18,17,18,13,8,3で合計101。
全非零D entryのweightがblock規則を満たすことも確認し、これらを独立な101方向の下界として使う。

同じ点で以下の全59列を再構成した。

- 各GL₅ elementary matrixに対応するEuler gaugeと実PGL微分をinterleaveした50列。
- scalar一列。
- 全七source微分。
- 新しい実curve方向R一列。

全150×59 entryは公開producer matrixと一致する。
49×49 minorは厳密に`239770014869548800000000≠0`。
全D×actual matrixは零である。
従ってker Dのdimensionは49以上、rank Dは101以上なので、両者は正確に49及び101。
大きなmatrixのgeneric rankを境界へ適用したものではない。
旧Fはrank48、Eulerは24、旧group/gauge/scale spaceは47も独立に再確認した。

## 3. 接空間の和をactual image下界と混同しない

新しい下界には**一つのjoint actual family**がある。
全λ・ε曲線へ、全六磁気normalizerの逆PGL変換とinverse tensor scalarを作用させる。
この同じ点でα=13/2、β=−15/2、I=1は全てunitなので、近傍でregularなjoint familyとなる。
ε=0では元の全七sourceそのもの。
したがってその七parameter微分とε微分は、上で保存した七方向及びRを同時に実現する。
更に実PGL、scalar、Euler gaugeのparameterを作用させる。
Euler gaugeはprojective Poisson構造の同じraw lift代表を変えるだけなので、全350 scheme内のactual方向である。

このregular familyの微分rankは49であり、標数零ではactual image closureのdimensionも49以上。
その像を含む点のlocal raw Poisson dimensionは49以上、一方全Zariski tangentが49なので、
local dimension=embedding dimension=49でregularとなる。
Euler gauge二十四とprojectivization一を除き、projective local dimensionは24。
regular local ringはdomainなので、この点を通る全空間成分は唯一である。

全quadratic_A sourceのMPS40閉和への帰属は既受理の別証明である。
この点を通る唯一の全空間成分は従ってMPS40の一つ、dimension24である。
特定の番号やgeneric orbit dimensionを追加推測しない。

## 4. 旧true48cutの支持原点目標は偽

旧weight-frameの真正法線N一つについてD N=0を再検査した。
rank[F48,N]=49と、Dの独立101標準座標列を合わせると、全150のdirect sumが得られる。
よって旧L=span(101 implicit標準列,N)は、F48に横断するdimension102の真の線形切片である。
これは新たに切片を選び直したものではない。

raw smooth germ dimension49に対し、旧48cutのdifferentialsは独立48。
従ってそのlocal cut supportはdimension一であり、dimension零ではない。
既受理のtorus固定IFTによる全18変数模型との同定を使えば、その旧模型の支持もdimension一。
そのnormalのpure powerをideal memberとして原点支持へ縮める旧目標は偽である。
旧Mora/unitの時間切れ結果を成功へ書き換えず、新しいactual curveと次元一致で数学的状態を訂正する。
旧idealの明示parameter表示又はscheme厚さは計算していない。

将来の別sourceの帰属へこのsmooth点を使うには、そのsourceの実族が実際にこの点を通ることが別に必要である。
退化因子型や未証明のλ橋の一致だけでは適用しない。

## 5. 固定証明書

finite job gate後0.518秒、最大RSS18,816KiBで全exact arithmeticが終了した。
producer及びCASはimportせず、既存の独立Fraction/QQ(λ) helpersだけを用いた。
新anchor、追加jet、局所ideal計算は無い。

|artifact|SHA256|
|---|---|
|code/classification_quadratic_a_smooth_review.py|e2f8c084f910e3792389a531869e46d3641dcce43da601b32a2164f6f02bb347|
|reports/CLASSIFICATION_QUADRATIC_A_SMOOTH_REVIEW_20261107.json|5de4cda347acbed15fe735dabafddf6d693fcffd8a79d778756cef1675106a79|
|reports/CLASSIFICATION_QUADRATIC_A_SMOOTH_REVIEW_JOB_20261107.json|a4af50d5dcfadd4941a7971185ad6fa06b6c296874583adad5ac907b30f40c34|

全元source、旧frame、literal曲線、既受理whole-source証明及び新producer packetのpath/hashはjobとreview reportに保存した。
旧artifact、global台帳、原稿及びGitは変更していない。
