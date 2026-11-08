# reduced triangle・m二double normal leading全枝：別作者の完全確認

2026-11-08、`/root/quartic_review`。入力は20261107凍結proof packet。
**全sourceの数学的還元、独立Fraction算術、全parameter補集合を受理する。**
nonconcurrent triangle R、rank二Poisson plane、actualD被約既約quintic、normal位数二、
generic double normal quadraticが対象。
h=f全源は空、h=L²M全源の像は23以下。
既独立leading排除と合わせ、全triangle normal-double m二G非優越が完了する。
h=L²Mの全点の既成分帰属、concurrent/nonreducedR/m三及び全G/Eは結論しない。

## 1. 不変入力と新算術

HF入口manifest f8830301…及びH21入口539c2bdf…を入力とした。
producer actualD還元は最新7fa13189…（raw36のu²追加mixed termsを含む版）だけを使用する。
旧版のraw支持を無断に使わない。
HF report f7dd10a2…、H21 report a78ea661…はそれぞれ実行producer全文を含む。
current mutable runner/jobsをinput identityとしない。

新checker `classification_triangle_m2_double_review.py` はproducer/CASをimportせず、
既独立Fraction sparse演算だけから全式を再構成した。
job `jobs/triangle_m2_double_separate_review.json` のgate PASS後、
**2.519秒、42,896KiB**で次を確認した。

* HF二つの全1075×150 Casimir/divergence matrix、全labelsと非零150minor、kernel零。
* H21全十値の元150P0、全350 Jacobi零、actualG2。
  全350×30/J1、350×9/J2、complete kernels、全行消滅、独立性、非零minor。
* Euler25＋Levi13＋scale一の全150×39、Euler rank24、Levi stabilizer三、filtered fiber九。
* 全bad-row map、全39actualD一次項、既独立universal quadratic表とのexact basis/全表一致。
  optional quadratic compatibilityは一切実行しない。
* F*とFgenの全gradient syzygy space AR0/AR1/AR2を、全四変数monomialsから再構成。
  kernelが指定T又はT,Zの全linear multiplesで尽くされること。
* normal orderとdegreeから許される全weight tupleを整数で列挙し、全十二rootから十壁を得ること。

成功report `CLASSIFICATION_TRIANGLE_M2_DOUBLE_SEPARATE_REVIEW_20261107.json`
に実行reviewer全文、全input hashes、各unitと判定を保存する。
算術は後述の全actualD還元を自動認証するものではなく、その還元も別作者として以下で確認した。

## 2. HF：全raw sourceからcanonical Casimir derivationへ

既corrected leading源はX=Y=(0,−y,z)、Z=diag(s,t,−s−t)、K=k u²、k=1+s。
unitはk−s一であってk自身でない。
全P1のbad u-rowはc v² partial_x、c非零。
c零なら全u-rowがu倍、actualDもuで割れるので指定既約源に入れない。
全raw36を許すとu-columnにはb_i u²項も入り得るが、これは下のspan(u,v)へ入るだけ。
全P2九も保持される。

row divergence R_i=Σ_j partial_j Pijは
Ru零、Rv=−2ku、Rx∈span(u,v)、Ry=−2y+span(u,v)、Rz=2z+span(u,v)。
base P_fはdivergence零、trace X/Y/Z零、P1 base-pairのbase derivativeはnormal constant倍、
mixed P1のnormal derivativeもnormal linear、P2 base derivativeは零だから、全sourceでこの式になる。
raw fullP自体をstrict Poissonと仮定しない。

row divergence(E wedge R)=−5RなのでPi=P+(E wedge R)/5がcanonical divergence-zero lift。
既canonical BV/Casimir恒等式からPi sharp dF零。
従ってdense u非零上でXi=(Pi sharp du)/uはFを殺す。
これはglobal Euler correctionをactual stabilizerと取り違えた主張ではない。

Xiのy/z部分はeigenvalues−3/5,+3/5にnormal(u,v)からのlinear off-blockを足したもの、
normal blockはu固定/v→−3ku/5のnilpotent block。
Sylvester equationsは固定nonzero eigenvalue差で常に解け、k零も含めactual base-from-normal shearsでoff-blockを除ける。
Xiを−5/3倍すると

`S+N, S=y partial_y−z partial_z, N(u)=0, N(v)=ku, N(x)=a v²/u+bv+du, a!=0`。

NはC[u,u^-1,v,x]上locally nilpotent、Sとcommute。
各nonzero S-weightの有限N-cyclic空間ではS+Nが可逆なので、F∈C[u,v,x,yz]。
k零ならnonzero coefficientを持つN=g partial_xからFx零となり、F2=unit*xyz u²に矛盾。
**k非零はここで初めて導かれる。** k零のsourceを黙って落としていない。

## 3. HF：全invariant ringと二actualDのactual inverse

k非零なら
`J0=u²x−(a/3k)v³−(b/2k)uv²−(d/k)u²v`。
actual v-shear、x-shear及びnonzero diagonal/cube-root choiceでJ=u²x−v³へ戻せる。
全変換はlinearでSplaneを保ち、逆shear/逆scaleが存在する。
rational v/u等をglobal PGL同値とは呼ばない。

T=yzとするとu非零上のN-invariantはC[u,u^-1,J,T]。
u最小negative powerのnumeratorをu零へ戻したとき、
(J,T)→(−v³,yz)はpolynomial ring間のinjective mapなのでcancellationは出来ない。
従って元のpolynomial ringとのintersectionは**C[u,J,T]**そのもの。

degree(u,J,T)=(1,3,2)より全degree五Fは
`Au5+B u²J+C u³T+D JT+E uT²` の五候補。
normal order一のuT²を除き、leading F2からD非零。
actual x→x+(C/D)uでu³Tを除き、nonzero scaleを合わせると
`F=J(yz+B0u²)+kappa u5`。
kappa零は可約。B0零と非零を分け、明示的なinverse diagonal/cube-root mapsで

`F_epsilon=(u²x−v³)(yz+epsilon u²)+u5, epsilon=0,1`

へ戻る。B0非零のmapはold x=(kappa/B0)X、old z=B0 Z、old v=rV、r³=kappa/B0、他恒等。
B0零はold z=kappa Z、他恒等。
両方とも全nonzero unitsを明記し、零枝を残している。
両Fはx-linearのcoefficientとconstant termのgcd一から複素数上で絶対既約でもある。

## 4. HF：二つの全150必要空間は零

canonical liftPiにはdivergence零とPi sharp dF零が必要。
全150斉次二次係数を未知数とし、Casimir五componentの全degree六monomials1050行とdivergence25行を使った。
全row/column labels、全1075×150係数を独立に再構成した。
両方でnonzero150minorがあり、rank最大150／complete kernel零。
従って非零canonical liftは存在せず、HF全rank四源は空。
追加Jacobi又はPfaffian条件を解く必要は無い。
このfull150空性は、二actualDへの**全source還元**と組み合わせて初めてwholeHF空性になる。

## 5. H21：canonical weightsと二fixed源の全空性

corrected leading源はDelta=diag(−1,2,−1)、Z=x partial_z、K=kappa uv。
row divergenceからRu=kappa u、Rv=−kappa v、Rbase=X+Y+span(u,v)。
5Xi/3のbase weightsはmu_j=s_j+4d_j/3+kappa/3、normal weights(u,v)=(0,−kappa)。
これをraw fullsourceから読むので、strict原代表又はgeneric off-block zeroは仮定しない。

fixed bad-yではweights(x,y,z,v;u)=(−2,4,−4,2;0)。
actual shearで全linear off-blockを消すと、resonant N=c v²/u partial_yだけが残る。
S,N commuteし、NF零が必要。
N(F2=unit*x²y u²)のnormal order三はnonzero*x²u v²で、higher項から相殺出来ない。
従って全源空。

fixed bad-zではweights=(−3,6,−2,−1;0)、同じ議論からFz零。
全weight-zero degree五monomialsはu5、u²x²y、xyv³だけ。
既約性とleadingから三係数全てnonzero、actual inverse diagonalsで
`F*=xy(xu²+v³)+u5`へ戻る。

F*のgradient syzygy degree≤2について、y零でb,cがy倍、次にx零でcの残linear部がx倍、
x=v零比較でそれが零になる。残式から全syzygiesは
`T=(3x,−6y,0,v)`のpolynomial multiples。
独立linear algebraでもAR0/1/2 dimensionsは0/1/4、全kernelがこの形と確認した。
vertex zを戻しても全quadratic Casimir bivectorのimageはTとvertex方向の二spanに入る。
従ってrank四は出ず、fixed bad-z全源も空。

## 6. H21：whole oneparameterを十壁へ絞る証明

bad-x源はY=diag(−2r−1,r+1,r)、kappa=3r+1、r!=−2。
canonical weightsは
`(mu_x,mu_y,mu_z,nu_v;nu_u)=(−r−2,2r+4,2r−1,−3r−1;0)`。
mu_i=0又はmu_i=nu_vの全linear resonanceは正確に{−2,−1,0,1/2}。
その外ではactual shearsでlinear off-blockを消せる。
bad termはmu_x−2nu_v=5rなのでr非零上でrational x'=x+lambda v²/u、lambda非零により除ける。
このrational shiftはinvariant必要形を読むためだけで、global PGL同値に用いない。

S=u x'=ux+lambda v²とする。
u最小negative exponentのnumeratorをu零へ制限してcancelするには、
固定y,z powersのS^a v^bとS^a' v^b'に2a+b=2a'+b'が必要。
weight零も使うと5r(a−a')零。r非零なので同じmonomial以外はcancel出来ない。
従って全Fは**C[u,S,y,z,v]**へ入る。

degree条件はq+2a+m+n+b=5、最小normal orderはq+a+b。
異なるmonomialsの最小normal項は一意なので、leading unit*yS²以外の全項はa+m+n≤2。
全nonnegative tupleを列挙すると、weight
`(−a+2m+2n−3b)r+(−2a+4m−n−b)`
の零は正確に十二値
`{−7,−2,−1,−3/4,−4/7,−1/2,−1/3,0,1/7,1/2,1,3}`。
checkerはこの全有限tupleとlinear resonance集合を整数/Fractionで再構成した。
monomial samples又はgeneric rankから推測したwallではない。

r零ではrational shiftを使わず、commuting nil部はc v²/u partial_x+d v partial_z。
F2=x²yu²のnormal order三への作用はnonzero*x y u v²で、z方向項では消えず全空。
r=−2はG2零となるm三以上の別義務。
r=−1,1/2のresonant linear chartは十壁へ戻し、off-blockをzeroにしたsliceで代替しない。

壁外ではF=unit*yS²+A u5、A非零。
actual scalesで `Fgen=y(ux+v²)²+u5` へ戻る。
Q=ux+v²のgcd argumentから全degree≤2 syzygiesは
`T=(2x,−4y,0,v)`、`Z=(2v,0,0,−u)`のpolynomial multiples。
独立linear algebraもAR0/1/2 dimensions0/2/8と全kernelを確認した。
vertex zを戻したimageは高々三次元、skew rankは偶数なので高々二。
従って壁外のwhole源も全空。

## 7. 十wholewallとactual group fiber

残十値では原P0と全30P1／全9P2へ戻した。
共通のfull39 quadratic表は既独立tangent receiptと全basis/全係数が一致し、再実行しない。
新しい全一次作用素とactualD一次項だけを各P0から独立再構成した。
complete kernelsと全minorにより次の全source上界を得た。

| r | k1 | k2 | 像上界 |
|---|---:|---:|---:|
|−7|7|0|22|
|−1|8|0|23|
|−3/4|7|0|22|
|−4/7|7|0|22|
|−1/2|7|1|23|
|−1/3|7|0|22|
|1/7|7|0|22|
|1/2|7|0|22|
|1|7|0|22|
|3|7|0|22|

Euler rank24。全150×39 actionから全Levi stabilizerは37−rankG三。
全六filtered shearsはassociated-graded P0を保ち、全higher源へ作用する。
そのgroup graphはPGL×sourceの各fiberにdimension三＋六を与える。
source上の固定点を理由にこのfiberを失わず、Leviとunipotentも重複しない。
固定P1上のP2は空又はkerJ2のtorsorの部分集合なので、全nonlinear sourceにも
`24−9+k1+k2` が有効。
全壁で23以下だからoptional quadratic compatibilityは不要で、**実行していない**。

## 8. 完了する義務と残余

HF全空、H21二fixed及び壁外全空、十wholewall像≤23。
既独立のmoving root／h=L³全P0排除と合わせ、
**reduced nonconcurrent triangle・m二double normal quadraticの全G像は23以下**。
有限constructible像の閉包も同じ次元上界を持つ。
H21の全pointwise E帰属、concurrent double/nonreducedR/m三/他surface及び全G/Eは別義務。
成分台帳、global register、Git又は主稿をreviewerは変更していない。
