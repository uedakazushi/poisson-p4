# R三重直線・m三squarefree全source空性の独立review

2026-11-07。producerとは別作者の数学的読解と標数零Fraction算術。**PASS**。
実退化因子が既約被約quintic、rank二restriction planeのzero cubicが3L、
plane沿い位数が三、generic normal binary cubicがsquarefreeである全sourceは空。
R=2L+M、m2/m4、m3の重根、他支持及び全G/Eの完了は主張しない。

## 1. 全二first connectionの幾何を独立に確認

producerの `PLANE_NONREDUCED_TRIPLE_LINE_M3_SQUAREFREE_FIRST_JETS.md` と
`PLANE_NONREDUCED_TRIPLE_LINE_M3_SQUAREFREE_EXCLUSION.md` を全読解した。
Hamiltonian associated-graded connectionは非被約R上の高位極を許し、log性を仮定しない。
U=P2\L=A2上のnonzero parallel cubicは一点でもmultiple rootになることが無いので、
root incidenceはfinite étale degree三。A2(C)のsimply-connected性から三つのdegree一coverへ分かれ、
normal baseのfinite birational mapとして各coverはA2自体である。
この議論でnormalization/basechangeや未証明のflatnessを仮定していない。

Pic(A2)=0より三rootをpolynomial vectorsで表せ、各pair determinantは nowhere-zero polynomial、
従ってnonzero constant。overall coefficientもdiscriminantからunitである。
generic cubicのbase degreeが高々二なので、三factorのbase degreesの和も高々二。
top homogeneous productはdomain上で零にならないため、見掛けのleading cancellationでこの次数を下げられない。
非constant factorが一つだけなら二つのconstant independent vectorsとのdeterminantsによりそれもconstantとなる。
従ってconstant三根、又はdegree pattern(0,1,1)だけが残る。

後者はconstant rootをzとして、他二rootのw coefficientsがconstant非零なので、
w+ell z と w+(ell+c)z、c≠0へ戻る。
ellのconstant項だけをnormal shearで消し、nonconstant affine-linear部分をactual base coordinate b/aにする。
base/normal transformationsはいずれも可逆で、三根相異性とD3零性を保つ。
twisted cubic actionのLie stabilizerは零なので、全parallel connectionが一意に復元される。
nonlog又はoff-diagonal接続のparameterを捨てていない。

結果の二raw P0は

\[
 P_{0,s}=3a^2\partial_b\wedge\partial_c,\qquad
 P_{0,u}=3a^2\partial_b\wedge\partial_c+3az\partial_c\wedge\partial_w.
\]

actual D2=0がnormal-normal quadraticの三係数を零にするので、未処理のraw Kも無い。

## 2. 全regular parallel spaceも別算術で確認

base-degree二・normal-degree三の全二十四monomialsを使い、scalarでは
a∂bG=a∂cG=0、unipotentでは(a∂b−z∂w)G=a∂cG=0を全48×24行列へ展開した。
完全核は四及び三次元で、順に

\[
 a^2\operatorname{Sym}^3(z,w),\qquad
 \langle a^2z^3,az^2(aw+bz),z(aw+bz)^2\rangle
\]

に等しい。全amplitudeを保持し、chosen cubicでsourceを切っていない。

## 3. 全30/9とactual leading imageの独立再構成

新checkerはproducer/producer helper/CASをimportせず、既に独立作成されたFraction sparse engineだけを再利用する。
weight一のraw36方向に対してcanonical30＋Euler6の全rank36を確認し、weight二raw9を全て残した。
producerのlast12のsource orderはcanonical orderの置換として明示した。
全350×30/9、全24×30 actualD3 map、全parallel-membership行をゼロから再構成し、
producerの全entryと一致を確認した。

| first connection | fullP1 constraint rank/kernel | fullP2 rank/kernel | 全actualD3像 |
|---|---|---|---|
| scalar | 6 / 24 | 0 / 9 | 0 |
| unipotent | 12 / 18 | 3 / 6 | span(a²z³) |

独自Fraction eliminationで全kernelを作り、producerの保存kernelの全row代入・独立性・nonzero rank minorも確認した。
全D3 mapを完全P1 kernelへ掛けた24×24又は24×18 image matrixを保存し、rank0/1を確認した。
unipotentの非零rowはa²z³だけで、independent rank-one minorは3。
各全P2と全P1 quadratic polarizationsのnormal-degree三coefficientsも零であり、
これは全P0+P1+P2のactualD3そのものである。

scalarではm三に必要なleading非零条件が失敗する。
unipotentでは非零leadingもtriple rootなのでsquarefree条件が失敗する。
任意のfull nonlinear Poisson solutionはこの完全linear必要源に入るため、全sourceで空性が従う。
高次Jacobiを解かずに得る最弱の支持結論であり、追加minor/GB/jetは行わなかった。
旧coarse reportのfalseと41/35上界は不変である。

## 4. 独立receipt

finite job `jobs/plane_nonreduced_triple_line_m3_review.json` のgate後に一回実行し、
9.931217秒、RSS17,840KiB、全独立算術PASS。
gate自体は数学の認証でない。

| artifact | SHA256 |
|---|---|
| code/classification_plane_triple_line_review.py | bbfda6fb7debbb3ab4e8b5d3bc5627547c4c208622aacb78ed32f25cb7c88755 |
| code/classification_plane_m4_special_review.py | 3c43fe85efcae39dbd0d2f8228271ce5fa291830315e62c9169ac5924b7408f8 |
| reports/CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_REVIEW_20261107.json | dc3ce248e18fcfffd8e9b411882bf87f257419ff765134bc9095ba2ec13b3384 |
| research/classification/jobs/plane_nonreduced_triple_line_m3_review.json | 37f64fc2bb61843d8540292e8a27c6fdfea0b3524982bfac5c4f326d53fb54da |

producerの新leading-exclusion receiptは計算終了後にhashを照合するだけで、解法として使わない。
全実行reviewer bytesと入力hashは独立reportへ保存した。
