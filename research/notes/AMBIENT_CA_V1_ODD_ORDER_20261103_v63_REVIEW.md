# CA V1のactual57制約とgeneric V全完成環の独立監査

2026-11-03。成功V2、actual odd rank9 supplement、および新しいreview-gated completion companionの適用鎖を承認する。
Hessian helperはAUTHOR contributionとして扱い、それ自体を独立reviewの根拠にしない。
CAのV1とその稠密V軌道が範囲であり、CB、B18原点、actual CA/CB交叉Dの全型は含めない。

## 1. 実sourceと独立の物理係数計算

先行する全proper CA旗模型、P0の全旗fiber、全critical schemeの独立receiptを固定した。
保存C20は実fixed-flag linear fiberである。全critical Schur150式へmonic V曲線(u=w=0)を代入し、
恒等零と両端P0,V1を確認した。全GL25作用からtarget orbit rank20を再計算した。
29source frameは実旗9・fiber20座標の中で可逆であり、group20と補fiber9は全normalizer sourceのetale座標になる。
GL列10は全150係数でP自身と一致する。これは非中央E22であり、中央GLによる作用ではない。

独立実装は、10次元外積index作用と15次元quadratic monomial反変作用を分けてGL Lie作用を構成する。
全J150×29、ordered exponentialの全contracted Hessian150×29を再構成し、author checkpoint/helperと一致した。
diagonalもL_i²Pを含む。checkpointに残るI+sEという旧firstjet recipeを高次で使わない。

## 2. actual57 systemと三次残余

拘束は20 orbit targetのeven値、29 quadratic-base targetのodd値、8 actual regular odd値である。
raw Y=0というsliceではない。全57 blockを独立に再構成して、保存inverseとの左右積がI57であることを確認した。
八oddとorbit20のregular minorは100/6561である。

独立のordered Lie exponential seriesで純kernel曲線と補正後のsource曲線を評価し、成功V2の全150有限seriesと比較した。
producerの5×5座標pushforward evaluatorは呼ばない。
actual odd8行は27,28,29,35,60,65,67,68で、29 base出力とこれら8出力はevenである。
他方、全残余の非零三次は正確にraw63=800/27 epsilon³、raw64=-40/27 epsilon³の二つである。
source kernel座標のodd部分はepsilonそのものに固定されている。

57 unit systemはepsilon²でrescaleしたpaired-source方程式の真の形式陰関数を与える。
二つのsource点は同じquadratic base値を持ち、全fiberの長さ2を後述のproper論法で確定するとconjugate pairそのものになる。
そのtracefree kernel座標はt=epsilon、nu=t²であり、非零三次はactual odd coefficientのnu一次classを検出する。
高い五次・七次を計算しなかったことは、この非零一次classの結論を妨げない。

## 3. 全actual odd firstjetとcritical scheme

Jrの28 target Y blockを逆にし、W=F_(kappa²)/2、Z=F_(kappa³)から、全actual odd firstjetを再構成した。
Y-elimination後のmixed termをE、branch quadraticをwc、branch unitをwとすると、
branch-even部分を除く操作M↦M-wc*M_branch/wを使う。
tracefree座標への変更はodd firstjetを変えず、raw branch座標での行列は
`(BY-Bnu*L | Bnu)`、L=Jr_branch*(Jr_Y)^(-1)となる。
この式は単なるraw-kernel cubicではなく、全target-base補正を含む。

全150×29行列は保存supplementに一致し、9次minorは-70000000/59049である。
全150へのfactorizationも確認した。全source critical cotangent rem(H)のrankも9である。

actual28 targetとquadratic branchからregular29底Rを取り、tracefree sourceをB=R[t]/(t²-nu)と書く。
実像はA=R directsum tIであり、Iは全actual odd coefficientで生成される。
targetがYとnuを含むので、source critical idealは(t,I)で、critical quotientはR/(I,nu)に等しい。
critical schemeは20軌道を含みcotangent dimension20なので、全局所環は正則20次元である。
従ってJ=(I,nu)はRの正則9生成idealである。
actual Iのlinear rank9はJ/mJの全classを張る。NakayamaによりI=Jとなり、Rの座標をI=(s1,...,s8,nu)に選べる。

## 4. 全proper fiberとgeneric Vへの論理

P0で全proper旗fiberが一台長2であることは先行独立receiptの全旗計算による。
properかつ全fiberが有限なので、P0を含むtarget開で写像はfiniteとなり、fiber長≤2の開を取れる。
これはGL不変な条件である。

実V曲線のGL飽和閉包は既約でP0を含む。P0の全critical次元20がその上界を与え、
V1の20次元orbitが下界を与えるので、この閉包は20次元でV1 orbitは稠密開である。
finite/長さ≤2開はこのorbitに会うためGL不変性でorbit全体を含む。
V1 source微分のcorank1により選択点のlocal fiber長は2以上で、全長≤2と合わせて全fiberは一台長2となる。
別の旗preimageや別branchを捨てていない。

ここで最初からglobal normalizationを仮定しない。
まず全finite smooth source Bとその実completed IMAGE A=R+tIを取り、上の非零odd Iを確認する。
c!=0をIから選べばt=(tc)/cがFrac(A)に入り、Frac(A)=Frac(B)となる。
従ってこの時点で初めて有限正則Bが実Aのnormalizationと結論できる。
これが新completion companionの非循環な適用順序である。

## 5. 全完成型とalgebraの帰結

従ってcone完成環は `C[[g20,s8,t²,s_i*t,t³]]`、normalizerは`C[[g20,s8,t]]`である。
自由Poisson scalarは上のGL orbitに含まれ、projective版はg19となる。
全kernelはb=(nu,s8)、Y_i=b_i*tに対する36 Koszul関係と45積関係、計81本で生成される。
積でR+sum RY_iへ還元し、R-basis(1,t)とregular列bの全第一Koszul syzygyを使うと逆包含まで従う。
全81のpullback零も独立に確認した。GB/eliminationは実行していない。

coneはdim29、edim38、depth21、multiplicity2、projective版はdim28、edim37、depth20である。
normalizer conductorは(s8,t²)B、target conductorは(s8,t²,s_i t,t³)Aで、B/A=C[[g20]]tである。
有限normalizationはt²,t³を含むelementary subintegral拡大で、seminormalizationもBとなる。
このgermのsingular/nonnormal/nonCM locusは滑らかなconductor supportに一致する。

これらをB18原点、mixed cubic全体、CB、実Dへ広げない。
初回timeoutや旧候補cusp(2,7)を成功V2の証拠として使わず、全proper fiberとactual odd idealの鎖によって結論する。
review code: `code/ambient_CA_V1_independent_review_20261103_v63.py`。
receipt: `reports/AMBIENT_CA_V1_ODD_ORDER_20261103_v63_REVIEW/`。
