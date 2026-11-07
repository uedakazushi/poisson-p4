# rank零plane・重複度二の全sourceは既知MPS又はG非優越

2026-11-07、主義務G。producer全source証明。後続二/三parameter算術の独立review前はglobal採録を待つ。

**命題。** actual退化因子が被約既約quinticであり、rank零plane Sに沿う一般重複度が二とする。
この全Poisson sourceは、既知MPS40成分の和と、ambient projective image dimension二十四以下の有限個の族の和に入る。
従って未知β優越成分の一般元がこの枝に属することは出来ない。
全点のMPS membership、新成分不存在、m三以上又は全Gの完了を主張しない。

## 1. 全first symbolの被覆

新modular補題`REDUCED_DEGENERACY_ZERO_SURFACE_ABELIAN_CONORMAL.md`と別作者reviewにより、
被約Dではconormal可換。したがってrank零planeの全quadratic sourceはweight零35、一30、二9である。
weight零のtraceless co-Higgs entriesはpairwise wedge零。
そのsaturated共通lineはO(0)又はO(1)で、導出と全Lax式は
`PLANE_RANK_ZERO_MULTIPLICITY_TWO_FIRST_JET.md`による。

traceless entriesがすべて零なら、Π0=T∧Euler_normal+K∂z∧∂wのimageは高々三次元で、
Π0²=0となる。leading normal quadratic Pfaffianが零なのでm二に反する。
従って本命題でscalar-only first symbolを新たに除外する仮定は必要ない。

O(1)のcoefficient map M:C³→sl2がrank三、rank二/M(p)零、rank二/M(p)非零の場合は、
それぞれ既存のfullrank/kernel/nonkernel whole-source証明と別作者reviewで全像二十三、二十二、二十三以下。
rank一ならM=ℓ(x)Aなのでtraceless fieldはglobal linear projective field V=ℓUとconstant sl2 matrix Aの積である。
よってO(0)枝と併せて全globalV supersourceへ入る。
Vのdivisorial zeroを含め、saturationを正しく分け直す理由でsourceを落とさない。

constant Aは非零semisimple又はnilpotentの二型。
全weight零Jacobiとactual inverseは二つのglobal-vector first-jet notesに保存されている。
normalssでは[T,V]=0、K=kzw。
normalnilではK=az²+bzw、[T,V]=−bV。
m二はdet(Euler_base,T,V)非零に等しい。

normalssのsimultaneous Cartan枝はactual C* degenerationがSNC log-symplectic pointへ達するので、
`PLANE_RANK_ZERO_CARTAN_SNC_RECOGNITION.md`及び別作者reviewにより全MPS40へ入る。
他の全abelian二planeはmixed、regular-nilcentralizer、common-kernelで尽くされる。
normalnil b非零は全Jordan/共鳴源が二十二以下。
b零は全commuting Cartan/mixed/regular-nil/common-kernelを保持する。
これらの低parameter二十一源は別作者受理済みのlow-parameter closureにより二十四以下。
二parameter四族はtwo-parameter/hyperbola/nil-line三certificateで全sourceが閉じた。
最後の三parameter二族を以下で処理する。

## 2. 全parameterで一様なsource次元

二つの残る全normal formsは

`SS: T=tH+sN, V=H+N, K=kzw, t≠s`、

`nil: T=diag(0,1,t)_tf, V=diag(0,1,s)_tf, K=kz², t≠s`。

H=diag(1,1,−2)、N=E01は係数行列、fieldはx1∂x0である。
両者ともjoint base centralizer Aはdimension二。
block Levi G0=(GL3×GL2)/Gmはdimension十二、tensor scalarを加えて十三。
全symbolを固定するsubgroupにはbase centralizer二、normal centralizer一、
relative block scalar kernel一が独立に含まれる。
normalssではdiagonal normal change、normalnilではw→w+h zがKも保つ。
したがって全parameterでorbit dimensionは九以下。

whole三parameter sourceはbase十二以下、whole surface二parameterは十一、
whole curve一parameterは十、fixed pointは九以下。
special periodでstabilizerが変わってもこの一様下界を使うので、normal-form parameterに条件を付けたことから
根拠なくimage次元を下げる議論ではない。
plane choice六とtensor projectivization一を戻し、full higher fiber upper dを加えると
ambient upperは順に17+d、16+d、15+d、14+d。

## 3. full三parameter sourceと全linear walls

`CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_THREEPARAMETER_20261107.json`は
両全QQ[t,s,k] sourceのfull350×30/9、actual shear、complete kernels/maximal minorsを保持する。
両generic sourceはweight一kernel六＝actual shear全六、weight二injective。
actual finite shearで全Π1を消すとweight二JacobiからΠ2も零となる。
generic D0はSSでbase linear因子のsquare、nilでz²因子を持つので、
被約既約quinticを持つGのgeneric sourceは空である。
これはformal truncation又は既知境界pointへのmembershipからの逆推論ではない。

SSの全minor壁は十四linear planesで、各whole QQ² operator、全rational curve、全exact finite rootsを処理した。
非線形surface又は未処理curveは残らず、全最大ambient upperは二十四。
nil Cartanでは九linear planesが同様に閉じ、全最大upperは二十三。
各plane内の全curve factorはlinear-in-one-variableであり、A(t)s+B(t)=0、gcd(A,B)=1の全graphとして処理した。
A零のcurve部分は空なので、分母を払っても不明な垂直成分を捨てていない。
三parameter reportは次節の二surfaceを未解決として`all_passed=false`を保持している。

## 4. nil Cartanの二surfaceはactual同一のwhole graph

残る二式は

`F1=st−3s+3kt+t`、

`F2=st+2s+3kt−3k−2t−1`。

actual base permutation x0↔x1とtensor scalar−1は

`(t,s,k)→(1−t,1−s,−k)`

を実現し、全Π0を正確に同じnormal formへ戻す。
全higher sourceへの可逆linear PGL作用も同時に作用するので、F1とF2の全sourceは同一の像を持つ。
full5×5 P0恒等式とF1変換=F2をliteralに検証した。

F1はt非零で

`k=−s/3+s/t−1/3`

となる。t零ではs零であり、全vertical k-lineがleading determinant零のm三以上へ入る。
F2の分母t−1のvertical complementもt=s=1であり同じである。
従ってm二の全surfaceはこの一つのrational graphで尽くされる。
全D1/D2/shearをt倍してrankを変えずQQ[t,s]へ戻し、全kernel/minorを処理した。

全九curve壁のうち、graph domain又はm二条件で全体が除外されるものもその理由を保存した。
残るwhole curvesはlinear又は

`s=t²`、`2st−s−t²=0`

という全rational graphsであり、すべてのexact rootと分母補集合を処理した。
generic surface上界は二十二、whole curves及び非除外pointsも二十二以下。
source domain外でcleared matricesが零となるt零の見かけの大きいkernelを、実surface上の点とは数えない。
証明書`CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_CARTAN_SURFACE_20261107.json`は全sourceでPASS。
これで前reportの二つのwhole frontierを閉じ、nil Cartan全三parameter sourceは二十三以下となる。

## 5. 結論・証明書・残余

全globalV m二sourceは、normalss Cartanの既知MPS40又はimage二十四以下で尽くされた。
§1の全rank/co-Higgs分岐へ戻して命題が従う。
これは未知β優越成分の一般元を排除する弱い結論であって、二十四以下の族の全membershipを必要としない。

後続exact証明書はTWOPARAMETER、HYPERBOLA、NIL_LINE、THREEPARAMETER、CARTAN_SURFACEの
`CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_*_20261107.json`。
各実行producer全bytesをimmutable保存し、reusable runnerの全旧FunctionDef ASTを毎回保持する。
旧二十三file manifestの可変九入力、whole二parameter完了時の可変十二入力もhash付きarchiveへ保存した。
全QQ算術/finite inverse/parameter complementの独立reviewは各scopeを区別して付ける必要がある。

m三/四の全scalar、pointwise-dependent T,V、common-image、leading-Pf零を復活させる全sourceは別義務。
normalnil b非零だけは別の全multiplicity bridgeがあるが、本命題から全高m完了を推論しない。
rank二restrictionを持つsingular/vertical surface、全E及び全分類の尽くしも本命題の射程外である。
