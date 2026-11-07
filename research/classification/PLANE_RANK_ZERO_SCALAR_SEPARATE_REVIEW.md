# rank零plane・traceless coHiggs零・T非零の全十五源：独立review

2026-11-07。主義務G。**全source幾何・全operator・全parameter壁 PASS**。
実既約被約quintic Dを持つrank四Poisson tensorがplane上零、traceless coHiggsが恒等零、
scalar trace field Tが非零なら、その全projective像dimensionは二十三以下である。
T0K非零の別証明は本命題の入力とせず、T0K0及び非零traceless dependent global-V源も残す。
全E帰属又は全G完了は結論しない。

## 1. 全scalar sourceはsl3×sl2、十一parameter

`REDUCED_DEGENERACY_ZERO_SURFACE_ABELIAN_CONORMAL.md` の全重複度補題も読解した。
被約Dではrank零hypersurfaceは存在せず、rank零surfaceの一般点でlocal modular flowはそのcomponentを保つ。
normal-normal一次bracketの二係数はlocal modular fieldのnormal成分そのものなので零となる。
ordinary double又はm二を要さず、planeでconormalは全域可換である。

planeを(z=w=0)とする。全raw quadratic係数でplane上零となるものは九十方向、
Euler冗長性十を除くと八十。conormal可換性はnormal-normalのnormal次数一の六方向を消す。
残るnormal weightsは零35、一30、二9である。
weight一raw三十六からEuler六を除くcanonical三十、weight二はraw九全て。
これは一部のfirst-jet型を選ぶための追加仮定ではない。

weight零35のmixed coHiggs三traceless entryを全て零とすると、残る**全**sourceは

`Π0=T∧En+q_S(z,w)∂z∧∂w`, `En=z∂z+w∂w`,

`T∈sl3`, `S∈sl2`, `q_S=det((z,w),S(z,w))`。
Tのmatrix trace零はglobal projective fieldの一意なliftの意味であり、coHiggsのtrace field T≠0と矛盾しない。
三十二coHiggs方向のうちscalar八とvertical K三を保つので、affine dimensionは十一。
q_SはS=[[s,t],[u,−s]]に対しuz²−2szw−tw²で、全binary quadraticとの一対一対応である。

全Π0はJacobi零、Π0²零。
これはleading D2が零、すなわちm≥3であることだけを示し、higher sourceの空性は意味しない。

## 2. actual作用・全逆・Jordan網羅性

base g、normal h、tensor scalar σは実block PGL変換で、正確に

`(T,S) -> (σ gTg^-1, σ hSh^-1)`

と作用する。normal determinantはq_Sとvertical bivectorの変換で相殺する。
独立checkerはgeneric二×二hのadjugateで分母を払った全q_S共変恒等式を検証した。
これは特定hのsampleではない。
全block変換はnormal gradingを保ち、全30/9未知係数を可逆に移す。
特定canonical basisを項別に固定する必要はなく、Euler quotientの全graded sourceを保存する。

T≠0の全Jordan型は次の五つ。
D=diag(1,1,−2)、N2=E01、N3=E01+E12とする。

|T型|affine stratum dimension|GL3 centralizer dimension|
|---|---:|---:|
|regular semisimple|8|3|
|λD、λ≠0|5|5|
|λD+N2、λ≠0|7|3|
|N3|6|3|
|N2|4|5|

trace零と三次Jordan分解でこの五つを尽くす。
表示次元は9−centralizer dimensionにeigenparameter数を戻したもの。
Sはsemisimple κH（κ≠0）、nonzero nilpotent E、zeroの三型で、affine dimensionは3,2,0。

S semisimpleならtensor1/κでHへ戻し、Tの全residual eigenparameterを残す。
S nilpotentならσ倍したEをnormal h=diag(σ^-1,1)でEへ戻せる。
平方根を使わず、全逆も同じLaurent式で定義される。
S0は任意の非零tensor scaleを保つ。

σN3はbase g=diag(σ^-2,σ^-1,1)、σN2はg=diag(σ^-1,1,1)で元へ戻る。
mixedのN2係数も後者で戻せ、diagonal D部分はそのgと可換。
これら全inverse-diagonal identitiesを独立Laurent算術で確認した。

従って固定十型はN3/N2とS0,E,Hの六組、D/D+N2とS0,Eの四組。
一parameterはλD/H、λD+N2/H（λ≠0）と、diag(1,t,−1−t)/S0,Eの二組。
後二組のdomainはt≠1,−2,−1/2。
三固有値相異なら少なくとも一つの固有値は非零で、それをfinite permutationで選んでtensor scale一へ直す。
t=0,−1を除外しない。

最後の二parameterはdiag(λ,µ,−λ−µ)/Hで、
δ=(λ−µ)(2λ+µ)(λ+2µ)≠0だけを課す。
λ0、µ0は含む。三つのδ0 lineはrepeated-semisimpleの一parameter sourceへ戻り、
原点T0は本定理の外に保つ。
λD+N2/Hのλ0は固定N2/Hへ戻る一方、λD/Hのλ0はT0である。
これでnormalizationの全零補集合を区別した。

## 3. 固定十型とN2/K0のwhole factor

全150 Π0、全350 Jacobi、actual Pfaffian零を全十型で再構成した。
全350×30 D1、350×9 D2のentryを独立Fraction算術で照合し、
complete kernelの全行代入、独立性、nonzero maximal minorを確認した。

|fixed T/S|ker D1|ker D2|restored linear image upper|
|---|---:|---:|---:|
|N3/0|8|3|22|
|N3/E|6|0|19|
|N3/H|6|0|20|
|N2/0|16|6|31、線形boundは未閉|
|N2/E|7|0|18|
|N2/H|6|0|18|
|D/0|6|0|16|
|D/E|6|0|18|
|D+N2/0|6|0|18|
|D+N2/E|6|0|20|

N2/0では、完全D1 kernelの全十六方向について、∂aを含まないpairは∂b∧∂cだけと確認した。
全P2はbase pairのみなので、全必要sourceは非線形Jacobiを課す前から

`Π=∂a∧V+F2∂b∧∂c`

となる。actualDは厳密に`F2*(w*Paz−z*Paw)`。
F2=0又はcubic factor零ならrank四でなく、両者非零ならquadratic×cubicなのでD既約性に反する。
従ってこの**whole Jordan sourceは空**。
固定reportの`all_ten_fixed_G24_closed=false`及び線形31は不変のまま保存し、数学的factorで閉じた。
残る固定九源は全像≤22。

## 4. 一parameter四源：全二十root certificate

QQ(tau)を別実装のFraction有理函数体で扱い、全matrix・complete kernel・minorを照合した。
selected minor productに明示的tauを掛け、零parameterも必ずroot packetへ保った。
全factorは一次であり、各全rootの全350行をexact QQへspecializeして新kernel/minorを検査した。

|whole source|affine symbol dimension|generic d / image|domain内wall最大image|
|---|---:|---|---:|
|λD/H|8|6 / 19|19|
|λD+N2/H|10|6 / 21|20|
|regularSS/S0|8|6 / 19|23|
|regularSS/E|10|6 / 21|20|

全四familyで各五root、合計二十certificateを照合した。
前二源のroot集合は0,1/3,2/3,−2/3,−4/3。
後二源は0,−1,1,−2,−1/2。
outside-domain arithmeticも保存して再検査し、その値を現在のJordan theoremへ混ぜていない。
S0のt=0,−1ではd=11だが、fixed ratioでsymbol base七となるため全像は23。
残る例外も同じ全source復元で閉じる。

## 5. 二parameter：全十二line・九whole line・五十七root

QQ(λ,µ)の全generic350×30/9 matrixから、独立sparse QQ[λ,µ] Bareiss法でminorを計算した。
各exact divisionのremainder零を確認し、数値interpolationは使わない。
minor product×δの全十二linear factorの積との恒等式を照合した。

`λ,µ,λ+µ,3λ+2,3λ+4,3µ+2,3µ+4,3λ+3µ−2,3λ+3µ−4,λ−µ,λ+2µ,2λ+µ`。

generic kernelはD1六、D2零、symbol base十一から像22。
δ0の全三lineはexcluded Jordan領域とexactに判定し、前掲の別sourceへ戻す。
残る全九lineのparameterizationがentire A1と同型であることを確認し、
その全operatorをQQ(tau)で再構成した。
各line genericはd=6、symbol base十から像21。
各line自身のminorの全五十七rootもcomplete kernel/nonzero minorで確認した。
二十一packetはδ0でoutside-domain、残る全pointはd=6、symbol base九から像20。

全有限intersection五十四組についても、δ0補集合又はlineのnonzero-minor patch又は明示rootへ入ることを確認した。
global generic、一line上のgeneric、その全rootを合わせて全regular C²を覆い、
交点や零固有値を無検査のgeneric点として捨てていない。
二parameter源全体の上界は22である。

## 6. 全source次元を戻す矢印

fixed plane上のsymbolは(T,S)そのものなので、affine base次元はdT+dS。
normalizationで使ったscale/eigenvalueをこの数へ既に戻した。
各symbol上の全higher tensorはP1三十・P2九を持つ。
weight一はD1 P1=0、固定P1上のweight二はD2 P2+Jac(P1)=0という同一linear operatorのinhomogeneous system。
そのfiberは空又はker D2のtranslateなので、全higher fiberはd≤ker D1+ker D2。
より高いJacobiはこのsupersourceを狭めるだけである。

plane choice六、tensor projectivization一を足し引きして、全像はdT+dS+5+d以下。
固定ratioは真正一条件、二parameterのlineは一条件、finite pointは二条件である。
正規化後のeigenvaluesは元の非零scaleとの比なので、原source上の対応するeigenvalue式もcodimensionを正しく持つ。
有限permutationやconjugation frameの選択は次元を増やさない。
normalizer groupの次元を別に引いたり、PGL24をさらに加えたりしていない。

固定十源≤22（N2/0は空）、一parameter全四源≤23、二parameter≤22。
よってT≠0のwhole scalar stratumは**全像≤23**。
T0K非零、T0K0、traceless非零のpointwise-dependent global-Vをこの結論へ混ぜず残す。

## 7. 固定証明書

finite job gate後80.952秒、最大RSS35,172KiBで全算術が終了した。
producer/CASをimportせず、既存独立Fraction・RF・R2 arithmeticを再利用した。
旧fixed→one→twoの全FunctionDef ASTも保存bytes同士で照合し、互換appendであることを確認した。
producer reportは一切変更していない。

|artifact|SHA256|
|---|---|
|code/classification_plane_scalar_review.py|fbae9fad765e31321a7379c5616ef2bb385aa3ca4eb6fd9c15c30e070964d87b|
|reports/CLASSIFICATION_PLANE_RANK_ZERO_SCALAR_REVIEW_20261107.json|af88ff6d5c3e309e4e31ecb2fc4bd397be854f1d1258bab1f4833b900242ab6d|
|reports/CLASSIFICATION_PLANE_RANK_ZERO_SCALAR_REVIEW_JOB_20261107.json|4177b13500db1f240efc96852d7aa7398d9e84da25280878ac9bce5c186f5adc|
|producer input manifest plane_rank_zero_scalar_manifest_20261107.json|dfdc26f620600ffccb2f33fa0e93fc175b2e080e95f18d270a9a7a00ce14b68c|

全exact入力・helper・producer四型artifactはjob/report manifestへ保存した。
前提の全重複度conormal補題は本文§1で別に読解し、三十五/三十/九sourceの全方向数を復元した。
global台帳・原稿・Gitは変更していない。
