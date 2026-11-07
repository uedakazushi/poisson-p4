# dependent globalV全multiplicity：Ore入口と全支持の独立review

**独立PASS。** 被約既約actual quinticを持つrank四projective tensorで、rank零planeの
traceless symbolが非零global field Vとconstant sl2 matrixの積、trace fieldがT=tVなら、
この全sourceのambient projective imageは二十三次元以下である。
全m二、全rank零plane、全G又は全点のMPS帰属を本命題だけから主張しない。
旧六source reportのfalse、三十四/三十五/三十三の粗上界は変更していない。

## 1. 六源と残三源の完全なOre表示

先行 `PLANE_DEPENDENT_SIX_SOURCE_SEPARATE_REVIEW.md` で、全二十一Jordan族から六への
必要係数・全例外・actual inverseを確認し、全六の350*30/9と全kernelを独立再構成した。
SS Hは上界二十、nil Hは二十一、nil N/T=N/K=z2は必要Cが全kernelで零となりG空。
残るSS N/T=tN/K零、nil N/T零/K零、nil N/T=N/K零をここで閉じる。

全P1 kernelと全P2九は、Euler gaugeにより
Pi=Q(y)+partial_a wedge W(y)、y=(b,c,z,w)、a非依存へ入る。
この表示をproducerの一つのsliceから推論せず、次の独立な完全basisで確認した。
f2(z,w)三係数は全てQのbc bracketに保持し、P2のab/ac六係数は全Wに保持する。
Wのraw二十係数の残Euler冗長性はW->W+ell(z,w)E_yの二方向。
W_zのz2、W_wのw2係数を零にすると、全W十八係数となる。
ellがb,cを含む場合はnormal weight零のbase-base項を生じるため、このpositive-weight冗長性ではない。

通常SSでは非共通QのB,E各二係数と共通W十二係数が全P1十六basisとなる。
全350 weight一membership、canonical Euler quotient三十座標への復元を独立確認した。
選択した二basis minorsは -(t+1)2 と -(t-1)2。
QQ[t]のgcd一により**全複素t**で十六本が独立であり、標本tのrankから全lineへ外挿していない。
元D1 kernelはt=+-1,+-3以外で十六次元、+-1/3でも十六次元なので全て一致する。
t=1及び3では新f cohomologyを含む十七basisのconstant nonzero minorを確認した。
actual z/w swapとa signの可逆involutionがt->-tを全P0で実現し、負critical sourceも全て復活する。
nil trace一では四Q cohomology+共通十二が十六、trace零では七+十二が十九で、
全D1 membershipとnonzero constant minorsが元完全kernelと一致する。
全P2九がnew f2三+W六にそのまま等しいことも確認した。

a次数を全projective Jacobiで比較すると、正確に
[Q,Q]=0、E_y wedge [Q,W]=0となる。
前者はa*partial_a wedge [Q,Q]の係数から得られ、後者の逆も同じ式で成立する。
global divergenceによるstabilizerの延長やnormalizationのflatnessは仮定しない。

## 2. strict Q Jacobiは全cohomologyの正確なidealである

全四変数strict Jacobiの係数を全て独立展開した。f2三係数も初めから残している。
通常SSの全係数はDelta=det(B,E)のQQ[t]倍で、quotient polynomialsのgcdが一。
従って全tでidealは(Delta)であり、3t+-1を反転していない。
nil trace一もidealは(det(D,E))。
nil trace零では全quadratic coefficient-rowspaceが正確に

    2det(D,B)+fE_b,   det(B,E)+fB_b,   det(D,E)-fD_b

の三次元spanに等しい。SS t=3は(Delta,fE_b)、t=1は(Delta-fE_b,fB_b)へ正確に一致する。
各式が全Jacobi係数から生成されることと、全Jacobi係数が指定式から生成されることを
rowspace equalityで確認した。選んだ必要方程式だけで全supportを仮定していない。

## 3. constant factorは全essential3 theoremへ入る

Q=v wedge Z、v非零constant、L_vQ=0ならPf(Q)=0でactual Dはa非依存。
Dは符号を除きE_y wedge W wedge Qの係数である。
v微分ではv(E_y)=vの項が消え、L_vQ=0なので、残るのはE_y wedge [v,W] wedge Q。
module式をvでwedgeするとこの項も零となる。
よってDはaとvの二つの独立constant directionsに依存しない。
既約非零quinticなら本質二以下はbinary factorizationで不可能なので、本質三である。
既受理essential3全source上界二十三を**globalのまま**使い、vやplaneの選択次元を重複加算しない。

通常SS、nil trace一のdet零は二vectorsのcollinearityからこの枝に入る。
nil trace零f零では三vectorsがcollinear。全vector零でもf2*db wedge dcをconstant factorにできる。
f非零でfirst column零なら全Qはpartial_c因子、かつc非依存となる。
SS criticalのf零/first-column零の各枝も同様。
Q零ならactual D零で直ちにG空なので、非零Dへessential3 theoremを誤適用しない。

## 4. 残る全normalizerとcentral extensionの零非零

nil trace零、f非零かつfirst column非零の場合、三strict式とPluecker identityから
4D_bB_b+E_b2=0。actual w->w+h*zは
D->D-hE-h2B、B固定、E->E+2hB、f固定を実現する。
B_b非零ではh=-E_b/(2B_b)でD_b=E_b零、c->c+s*bでB_c零。
残る三式からD零、E=(0,-f)となりB型へ入る。
B_b零ではconicからE_b零、D_b非零。c-shearでD_c零とするとB零、E=(0,f)となりA型へ入る。
両first-entryとfの非零scaleは(a,b)のpaired scale、c-scale、tensor scalarで実現し、全inverseが可逆。

SS t=3ではf非零、B_b非零以外はconstant-factor枝。
残る式はE零を強制し、c-shear/scalesでB=(1,0),f一のS3へ入る。
SS t=1ではf非零、E_b非零以外はconstant-factor枝。
B=(0,-f)が強制され、同じactual changesでE=(1,0),f一のS1へ入る。
負critical periodsは前述のactual involutionで覆う。

全central f2も保持する。actual b->b+u(z,w),c->c+v(z,w)により
f2->f2-h_b*u-h_c*v+X_b(v)-X_c(u)。各canonical型でのimageは順に
A:span(z2,zw)、B:span(z2,w2)、S3:span(zw,w2)、S1:span(z2,w2)。
残る一係数gammaは、actual(a,b,c)->(lambda^-2*a,lambda^-1*b,lambda^-1*c)と
tensor scalar lambdaでgamma/lambdaへ移る。
これはP0とanchorを保ち、全W sourceを全W sourceへ戻す可逆変換である。
gamma零と非零を別に保持するため、全supportは四anchor×delta0/1の八whole Q orbitへ入る。
exact orbit dimensionを必要とせず、全W moduleに入る包含だけで以下の空性が従う。

## 5. 八つの全affine moduleはrank四既約Dを持たない

全W十八変数のfull350 projective Jacobiを独立に構成した。
strict Q Jacobi零から、a指数零かつquadにaを含むtrue140行だけが残り、他全行は恒等零。
非零minor、full kernelと独立free rows、全350 particular residual又はaugmented nonzero minorを保存した。
共通partial_aの方向同士の外積は零なのでactual DはWに厳密にaffine linear。
particularと全kernel方向のD係数を直接展開し、共通factorを確認した。

| Q | delta | rank | affine k | 補助像上界 | 全D factor又は不整合 |
|---|---:|---:|---:|---:|---|
| A |0|11|7|22|z3 |
| A |1|12|6|22|z |
| B |0|13|5|21|z2 |
| B |1|13|空|空|augmented rank十四 |
| S3 |0|11|7|23|w3 |
| S3 |1|12|空|空|augmented rank十三 |
| S1 |0|14|4|20|zw2 |
| S1 |1|14|4|21|zw2 |

従って八whole族は全てG空。factorは全affine solutionを割り、特定Rやhigher点だけの因数分解ではない。
それをactual inverse正常化で戻しても、proper factor又は不整合性は保存される。
constant-factor枝はglobal essential3上界二十三、最初のH二源は二十/二十一。
これで全dependent globalV枝の二十三上界が従う。

## 6. 固定receiptとproducer公開後の照合

初期六源の独立reportは `CLASSIFICATION_PLANE_DEPENDENT_SIX_SOURCE_REVIEW_20261107.json`。
Ore全入口は `CLASSIFICATION_PLANE_ORE_ENTRANCE_REVIEW_20261107.json`、0.273秒、
checker `classification_plane_ore_entrance_review.py`。
八moduleは `CLASSIFICATION_PLANE_ORE_MODULE_REVIEW_20261107.json`、0.650秒、
checker `classification_plane_ore_module_review.py`。
各jobのgate後にFraction/sparse算術だけを使い、producer/CASをimportしていない。
初期の任意二minorは共通(t+1)2を持ったため、actualswapに対応する+-1の二minorへ交換して再gateし、
最終実行で全basis/全strict idealを再確認した。失敗した選択を全t証明とは呼んでいない。

この二算術はproducer最終batchより先行した。
公開後、immutable `CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_ORE_20261107.json`、SHA256
`ce5bb0eefeebd47bbd6ecd75d0c1766e8c911c20e477bde7d143f1afabcc8ad2`、
executed producer bytes `e7df297669b1e9de3118d3a9d7f3f1b21a6325da31d0dfb5e0ab018961e6bad7`へ固定した。
全八canonical Q/W0のliteral一致、全350*18 affine matrixとRHSの一致を検査した。
producer列から独立列へのpermutationは
[0,1,2,3,8,9,10,4,5,6,7,11,12,13,14,15,16,17]、Schouten conventionのscalarは二。
全rank、k、不整合、像上界の表も一致した。
producerの特定のforward/back変換行列を信じる代わりに、独立な完全basisと全t minorsで
同じ全source等式を証明した。producer全中間liftを逐一再実行したとの主張ではない。
mutable runner/current old jobsには依存しない。

最終hashと全依存は `inputs/plane_dependent_ore_review_manifest_20261107.json`。
旧false、全111/curve/ambient packet、初期六source及び今回の独立reportは凍結のまま保持する。
