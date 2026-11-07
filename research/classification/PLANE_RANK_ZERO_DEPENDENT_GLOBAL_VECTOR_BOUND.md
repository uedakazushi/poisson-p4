# rank零plane・dependent global field全支持のG非優越

2026-11-07、主義務G。producerの全source証明と有限QQ証明書。
旧six-source reportの`all_passed=false`と粗い上界は変更しない。
全六kernel及び八fixed moduleの数値は、producerをimportしない別作者の独立Fraction計算と一致する。

**命題。** 被約既約actual quintic Dを持つrank四projective Poisson tensorを考える。
rank零planeに沿う可換conormal first symbolのtraceless部分が、非零global projective field Vと
constant sl2 matrixの積であり、trace field TがVに比例するとする。
この全sourceのambient projective image dimensionは二十三以下である。
重複度二を仮定しない。従って未知β優越成分の一般元はこの枝に属さない。
全点のMPS membership又は全rank零plane/Gの尽くしは主張しない。

## 1. 全Jordan型から六sourceへ戻る

`PLANE_RANK_ZERO_DEPENDENT_GLOBAL_VECTOR_FIRST_JETS.md`は、Vの全五Jordan型、
regular semisimpleの固有比、Tの全比例係数及びnormal binary bracketを保持する。
`PLANE_RANK_ZERO_DEPENDENT_GLOBAL_VECTOR_HYPERPLANE_REDUCTION.md`は、
full Π1の極端normal monomialに対するprojective Jacobiを使い、全支持を六sourceへ還元した。
その逆normalizationsはactual base/normal linear changesとtensor scalarである。

normalssでC又はDvが零、normalnilでCが零なら、実hyperplaneがPoissonとなりDはlinear factorを持つ。
この補集合をG空としてから、全係数式がVをdouble semisimple H又はrank一nilpotent Nに制限する。
normalssではHの全sourceはT零/K=−5zw、Nの全sourceはT=tN/K零となる。
normalnilではH/T零/K=(5/2)z²、N/T零/K零、N/T=N/K=z²又は零の四sourceを保持する。
trace非零の場合にK非零を誤って落とさない。

旧`CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_DEPENDENT_20261107.json`の全350×30/9から、
H二sourceの全像はそれぞれ二十、二十一以下である。
N/T=N/K=z²のfull weight一kernelはCを強制的に零とするのでG空。
この三sourceについて、linear必要条件だけで得た十分な結論をそのまま使う。
残るN/K零の三sourceを以下で閉じる。

## 2. 全sourceをa非依存Ore形へ移す

座標を(a,b,c,z,w)、N=b∂aとする。
`PLANE_RANK_ZERO_RANK_ONE_NILPOTENT_ORE_REDUCTION.md`の§1–2によれば、
全weight一kernelとfull weight二九係数をEuler gaugeで

`Π=Q(y)+∂a∧W(y), y=(b,c,z,w)`

へ移せる。Qのstrict Jacobi及びWの全affine module条件は正確に

`[Q,Q]=0, E_y∧[Q,W]=0`。

新reportの`full_source_Ore_coverage`は、通常SS、SS t=1、SS t=3、nil trace一、nil trace零の
五full kernelsからこの変換をliteralに再構成し、各全source係数間の変換行列の行列式が**−1**であることを確認する。
従って新しいQ/W表示はchosen sliceではなく、全Π2を含む全sourceの可逆表示である。
WのEuler二冗長性を固定した全十八係数も保持される。
これはglobal divergenceに由来するstabilizer延長の仮定を使わない。

Qがnonzero constant field vをfactorに持ち、L_vQ零なら、PfQ零でactual Dはa非依存である。
module条件をvでwedgeするとv(D)零も従う。
従ってDは独立な二constant directionsに依存せず、本質高々三変数となる。
既約quinticなら本質二以下はbinary factorizationに反するので本質三。
別作者受理済み`ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md`の全source上界二十三を適用する。
この議論は大きなrank二linear subfamilyの次元を数えることを必要としない。

## 3. strict Q Jacobiの全支持はconstant factor又は八canonical Q

Ore還元note§3–6の全normalizerとcentral extensionを使う。
通常SSでは二vectors B,EのdeterminantΔに対して
`(3t+1)Δ=(3t−1)Δ=0`が必要なので、全tでΔ零。
nil trace一でも全二vector determinantが零となり、いずれも前節のconstant factor枝に入る。
SSの通常full kernelの例外t=±1,±3は全て復活させ、actual normal permutationで負periodを正periodへ移す。
t=±1/3はminorの見かけの根であり、元の全kernelが通常basisに等しいことを旧reportが確認している。

nil trace零の全三vector/f source、SS t=1及び3の追加一cohomology及び全binary central extension f2を保持する。
strict Q Jacobiの正確なquadratic ideals、nil conicの両chart、f零/first column零の全補集合、
actual normal/base shears、central coboundaryの全image、残るcentral amplitudeの零/非零を使うと、
constant factor以外の全支持はA/B/S3/S1の各δ=0/1という八Qへ移る。

新reportの`strict_anchor_ideal_certificates`は全80 strict Q Jacobi係数と指定idealの双方向liftを保存する。
SSではQQ[t]のBezout恒等式を保存し、divisionで見えないperiodを捨てない。
nil Plücker→conicのlift、normal/base shear、全central coboundary、γ→γ/λのactual diagonal変換も
全parameter恒等式として検証している。anchorを選んだ八点から全支持を推論する議論ではない。

## 4. 全W十八変数のexact affine moduleとactual D

各canonical Qに対して、元の全350 projective Jacobiとmoduleの全140式は同じaffine rowspaceを持つ。
全十八W変数へのmatrix、constant term、非零minor、full kernel、particular solution又はcokernel contradictionを保存した。
compatibleな場合、actual DはWにaffine linearなので、そのparticular及び全kernel方向のD係数の共通因子が
全affine familyのactual Dを割る。必要式を満たした部分だけに限定した因数分解ではない。

| Q型 | δ | matrix rank | augmented rank | 全W fiber dimension | 全actual Dの共通因子 | 結論 |
|---|---:|---:|---:|---:|---|---|
| A | 0 | 11 | 11 | 7 | z³ | G空 |
| A | 1 | 12 | 12 | 6 | z | G空 |
| B | 0 | 13 | 13 | 5 | z² | G空 |
| B | 1 | 13 | 14 | — | — | 全Jacobi不整合 |
| S3 | 0 | 11 | 11 | 7 | w³ | G空 |
| S3 | 1 | 12 | 13 | — | — | 全Jacobi不整合 |
| S1 | 0 | 14 | 14 | 4 | zw² | G空 |
| S1 | 1 | 14 | 14 | 4 | zw² | G空 |

不整合二族ではexplicit cokernel combinationのconstantがそれぞれ−4、−24となる。
他六族は全affine solutionでD零又は上のproper linear factorを持つため、被約既約quintic Gには入らない。
actual normalizationsを逆に戻してもproper factorと不整合性は保存される。
補助的に全orbit/source/plane/scalarを戻したimage上界も22,22,21,23,20,21で十分だが、
主結論にはより弱い因数分解だけを使う。

## 5. 結論・証明状態・残余

六sourceのうち三つは§1、残三つのconstant factor部分は§2の本質三全source上界、
それ以外は§4の全G空で閉じる。よって命題の二十三上界が従う。
full350だけのscheme-theoretic分類、局所唯一性、全成分帰属又は新成分不存在は必要としない。

新有限証明書は`CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_ORE_20261107.json`、
jobは`jobs/plane_rank_zero_global_ore.json`、reusable runnerのmodeは`--global-ore`。
7.272秒、最大78,092 KiBで全検証PASS。
producer全bytesをreportに埋め込み、以前の全25 top-level FunctionDef ASTを不変に保った。
直前deb217f source及びcurrent dependent jobは`inputs/plane_rank_zero_dependent_frozen_inputs_20261107.json`から復元出来る。
全六sourceと八moduleの別作者算術reportは`CLASSIFICATION_PLANE_DEPENDENT_SIX_SOURCE_REVIEW_20261107.json`及び
`CLASSIFICATION_PLANE_ORE_MODULE_REVIEW_20261107.json`であり、本表と全て一致する。

Π0零のrank零plane、全scalar trace源、全Gの他のrank二/特異支持は本命題の射程外。
それらが別packetで閉じていても、本命題だけから全plane又は全G完了を推論しない。
旧false receipt、全m二certificate、全mutable-input archives及び独立review manifestのbytesは変更していない。
