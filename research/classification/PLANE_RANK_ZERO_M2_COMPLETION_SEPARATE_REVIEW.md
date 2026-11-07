# rank零plane・multiplicity二の全global-vector枝：別作者review

**独立PASS。** 対象は主義務Gであり、Eの全点帰属を主張しない。producerの
最終29-file manifestと13-job archiveを入力とし、以下の全source幾何と
独立exact算術により命題を受理する。

## 1. 命題と全sourceへの戻し方

実退化因子Dが被約既約quintic、Sがrank零planeで、Sの一般点における
Dの位数が二の場合だけを扱う。既受理のabelian conormal補題により、
全projective sourceはweight零35、一30、二9を保持する。
traceless co-Higgsの三entryはpairwise wedge零なので、共通imageの
saturationはO(k)である。P2上のsaturated rank一subsheafはreflexiveであり
line bundleとなる。非零sectionとH0(TP2(-2))=0からk=0又は1。
後者のcoefficient mapのrank三、rank二/kernel、rank二/nonkernelは既存の
別作者reviewを参照し、今回再実行しない。

rank一ではM=ell*A、したがってV=ell*Uはglobal projective linear fieldである。
Vのdivisorial zeroも許すglobal-V supersourceへ入り、飽和の取り直しで
sourceを落とさない。traceless部分零又はT,Vがpointwise dependentなら
leading quadratic Pfaffianが零なので位数二には属さない。

global-Vのconstant sl2 matrixはsemisimple又はnilpotent。
前者の全weight零Jacobiは[T,V]=0、K=kzw、後者は
K=az2+bzw、field bracket[T,V]=-bVである。
matrix unit Eijに対応するfieldはxj*d/dxiであり、matrix commutatorと
field bracketの符号を区別した。trace-free liftへ戻す際のEuler補正はKへ
吸収でき、Kは各表で全係数を保持する。

nilpotent normalかつb非零の全Jordan/共鳴族、低parameter21源、
semisimple normalのCartan全period帰属は既受理の証明を参照する。
新算術の対象は四つの二parameter源と二つの三parameter源である。

## 2. normal formの完全性と一様な次元復元

A=span(T,V)がregular matrixを含むなら、sl3 centralizer二次元から
Cartan、mixed、regular-nilcentralizerの三型。regular matrixを含まない場合、
非nilpotentな非regular matrixはdouble-eigenvalue semisimpleであり、
独立な可換matrixとの一般和がregularとなるので不可能。
残る非零matrixはすべてrank一nilpotentである。
rank一matrixの二次元線形空間はcommon image又はcommon kernel。
common imageはpointwise proportionalで位数二に反し、common kernelだけが残る。

H=diag(1,1,-2)、N=E01、J=E01+E12とする。
mixed normalizerはHを固定しNを任意の非零scalarでscaleする。
regular-nilcentralizerのactual GL3 normalizerは
J -> lambda*J+mu*J2、J2 -> lambda2*J2を実現する。
具体的にはJをsuperdiagonal一としたとき、上三角行列
G=[[lambda2,mu,0],[0,lambda,0],[0,0,1]]がGJ=(lambda*J+mu*J2)Gを満たし、
detG=lambda3非零である。逆polynomialは
lambda^-1*J-mu*lambda^-3*J2である。
producer両first-jet noteの非零係数による全normalizationを読解し、
零係数の各補集合が別表の行へ戻ることを確認した。
nilpotent normalのb零ではT,Vの独立scaleを使えるが、Vをnilpotentと仮定しない。
Cartanの場合は両fieldに共通する非零eigenvalue差の座標対を選べる。
三つの対のうち各非constant fieldで零となる対は高々一つなのである。

残る三parameter源は

* SS: T=tH+sN、V=H+N、K=kzw、t!=s。
* nil: T=diag(0,1,t)_tf、V=diag(0,1,s)_tf、K=kz2、t!=s。

全parameterでjoint base centralizerは少なくとも二次元、normal centralizerは
一次元、relative block scalar kernelは一次元で、これらは独立である。
後者のnormal centralizerはSSではnormal diagonal、nilではw -> w+h*zであり、
Kも固定する。Levi(GL3*GL2)/Gmとtensor scalarの全次元十三から
各normalized first symbolのorbit次元は一様に九以下。
したがってwhole三parameter、surface、curve、pointのaffine symbol base上界は
それぞれ十二、十一、十、九。plane choice六とprojectivization一を戻すと
全higher fiber上界dに加える数は十七、十六、十五、十四である。
これはparameter壁の全orbitに対して同じstabilizer下界を使う議論であり、
壁が群不変であるとの仮定を必要としない。
ここで引いている四次元はfirst-symbolのstabilizerであって、full tensorの
stabilizerと仮定していない。全sourceからfirst-symbol空間への射影を使い、
固定symbolの全higher fiberは群作用で同型となるため、base次元とfiber上界を加える。

## 3. 独立算術の最弱certificate

checkerはproducer/helperをimportしない。既存別作者のFraction/sparse-polynomial
関数だけを再利用し、各P0を上のT,V,Kから組み直す。
全150係数、全350 projective Jacobi、実Euler-signed Pfaffian、
全350*30/9 operator、actual shear30*6を再構成して保存値と照合する。

一般のG上界に必要なのはrankの下界である。各保存minorを独立にexact計算し、
非零と式の一致を確認する。選択minorの積の全factor identityを確認し、
全linear plane、全linear/rational curve、その全有限例外根を保持する。
全kernelの正確性を必要としない箇所でexact rankやschemeを余分に主張しない。
weight一kernel上界とweight二kernel上界の和はfull nonlinear fiberの上界であり、
weight三以上のJacobiを捨ててもこの必要条件による上界は有効である。

curve A(t)*s+B(t)=0のgraphではgcd(A,B)=1をQQ[t]で確認する。
A=0にcurve点はなく、分母を払ったmatrixの見かけの点を本来のcurve点と混同しない。
parameterizationが一座標をtそのものに保つことも確認する。
全factor productと全一次factor根の照合により、有限標本の検査にはしない。

## 4. SS regular-nil K零lineの真正余核

このwhole lineだけは粗いlinear boundが足りないため、別に全sourceを扱う。
complete kernel、actual shear、全first-weight slice、basis minorを照合する。
全slice二次Jacobi350行を直接展開し、D2の実pivot inverseから
全350の真正cokernel compatibilityを再構成する。

一般点ではshear rank六、slice二、D2kernel三。
選択二次式は4u0^2/(t+1)^2と12u1^2/((t-1)(3t-1))で、
binary resultantは2304/((t-1)^2(t+1)^4(3t-1)^2)非零。
rank、basis、inverse、resultantの全零点/極はt=+-1,+-3,+-1/3へ戻す。
各点の選択式は最初がmonomialであり、全選択式の共通variable因子が零個であることから
gcd一、従ってcompatibility idealのheight二以上を証明する。

t=+-1ではshear rankは**五**、slice四である。
全shear六を無料に引かず、actual additive shear groupのrank五の線形補空間だけで
P1をsliceへ戻す。残るshear kernelのhigher作用は正規化sourceに保持する。
一般上界二十四、t=+-3,+-1で二十四、t=+-1/3で二十三となる。
全P2kernel三を保持し、特定higher係数を選んだ点だけの検査ではない。

## 5. hyperbolaと二Cartan surface

nil mixedのcurve eta*tau-tau+3=0はeta=1-3/tau、tau!=0。
tau零の補集合は空、tau一はT,V依存で位数三以上である。
全operator/shearをtau倍したrank同値なQQ[tau]行列と全根を照合する。

nil Cartanの二surfaceは
F1=st-3s+3kt+t、F2=st+2s+3kt-3k-2t-1。
actual base swap x0<->x1とtensor -1が全P0で
(t,s,k)->(1-t,1-s,-k)を実現し、F1をF2へ送る。
これは全higher sourceにも作用する可逆PGL変換なので、一方の全sourceで十分。
F1はt!=0でk=-s/3+s/t-1/3。
t零ではs零、F2の対応するt一ではs一であり、全vertical k-lineが位数二の外である。
一つのwhole graphの全QQ2行列、九curve、全点を調べる。
s=t2、2st-s-t2=0もwhole graphとして扱う。

## 6. 固定入力・残余

producer最終source60e84e753c93f3f43a25ee9c884514552e8adc07423a4009eae22183643b818bは
immutable CARTAN_SURFACE reportのproducer_bytesから読む。
後続で可変となるrunner及び十三jobの現在bytesには依存せず、
manifestの当時hashとfrozen_jobs archiveの全十三bytesを照合する。
旧TWOPARAMETER/THREEPARAMETERのfalseは変更せず、それぞれの未閉line/curve/surfaceを
後続receiptで閉じる。低parameter二十一源の算術は既受理であり再実行しない。

位数三以上、pointwise-dependent T,V/common-imageを含むleading-Pf零の全source、
surface上のrank二restriction、全E及び全Gの残りはこの命題では閉じない。

## 実行receipt

独立checker `code/classification_plane_m2_completion_review.py`、SHA256
`34bf58aabb02f40045bf8fbb4fa08849efb9b1d3377a6794464b315070c69ef1`。
有限jobは `reports/CLASSIFICATION_PLANE_RANK_ZERO_M2_COMPLETION_REVIEW_JOB_20261107.json`、
SHA256 `906c1d1c6a82bbc656a49d1afddb25b0e83abf1f8b6e01d1fdbd1652cbeb3e1b`。
gate PASS後の最終全枝再実行は84.533129秒、最大RSS204028KiB。
結果 `reports/CLASSIFICATION_PLANE_RANK_ZERO_M2_COMPLETION_REVIEW_20261107.json`、
SHA256 `13312c0b204347262cb2a43c5f91dd454fead298d74b5ad5530deffdfcd41c2d`。

四QQ2源、二QQ3源、全二十三plane、二百三十二curve壁の階層被覆を確認した。
そこに記載された有限root955件と、SSnil K零lineの別の六root receiptを確認した。
これらはchart間の重複を含む出現数であって、互いに異なる点の個数ではない。
hyperbolaとnil K零line自体は上のcurve壁階層の後続receiptとして別に保持する。
SS MIXのwhole上界最大二十四、nil Cartanで二十三、Cartan二surfaceのgraphで二十二。
domain外のcleared行列に由来する大きい粗上界と、nil K零lineへ委譲した粗い失敗上界も
結果内に保持したが、それを実m二sourceの上界として使用しない。

最初の未完verifier実行は、同じ零・定数係数の特殊化でQQ(t)のgcdを何百万回も
繰り返す律速を除くため、証明書出力前に停止した。
Fraction Horner/cached特殊化とconstant Fraction determinantへ置換して再gateし、
上記最終実行で全枝を最初から確認した。producerの式、数学的scope、旧証明書は変更していない。

全m二への既受理依存は `PLANE_RANK_ZERO_MULTIPLICITY_TWO_FIRST_JET.md`、
`REDUCED_DEGENERACY_ZERO_SURFACE_ABELIAN_CONORMAL.md`、
`PLANE_RANK_ZERO_FULLRANK_REVIEW.md`、
`PLANE_RANK_ZERO_RANK_TWO_KERNEL_SEPARATE_REVIEW.md`、
`PLANE_RANK_ZERO_NONKERNEL_SEPARATE_REVIEW.md`、
`PLANE_RANK_ZERO_CARTAN_SNC_REVIEW.md`、
`PLANE_GLOBAL_LOWPARAMETER_CLOSURE_SEPARATE_REVIEW.md`。
これらのmath/source算術を再計算したとは主張せず、今回の全source橋へ適用した。
最終stable dependency manifestは
`research/classification/inputs/plane_rank_zero_m2_completion_review_manifest_20261107.json`。
