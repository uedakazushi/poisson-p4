# 既知142成分と七つの既約閉像による全Poisson支持の被覆監査

2026-11-08、PRIMARY E、独立合成担当 `/root/log_manuscript`。
入力は採録済みclaim1135まで、成分台帳142（一般rank四136、rank二6）で固定する。
主原稿・全点frontier・台帳・旧producer/review/certificateは変更しない。
新しい算術、operator、point、minor、profile、jet又はCASは実行しない。
これは全source定理の**点別の結論**を接続する純証明であり、Gの像次元上界からEを推測するものではない。

## 1. 結論と既知成分の入力

複素数体上の射影四空間の射影Poissonスキームの被約台を \(\mathcal P_{\rm red}\) とする。
採録済み142個の相異なる全空間既約成分の閉和を \(K\) とする。
各成分の成分性・相違は既証明を入力とし、今回再証明しない。
142個の同一性と依存claimの射程は添付projectionで凍結する。

以下の七つの既約閉集合を定義する。閉包・作用は全て \(\mathcal P_{\rm red}\) 内で取る。

| 記号 | 全actual source又は軌道 | 射影次元上界 |
|---|---|---:|
| \(C_G\) | Ga9 linear-wのcanonical軌道の閉包 | 23 |
| \(C_N\) | `node_flex_tangent_through_flex` の全六磁気源のPGL5飽和閉包 | 23 |
| \(C_F\) | `cusp_flex_secant_through_flex` の全六磁気源のPGL5飽和閉包 | 23 |
| \(C_C\) | `cusp_cusp_tangent_cusp_secant` の全六磁気源のPGL5飽和閉包 | 23 |
| \(C_D\) | constant-modular最高double-lineのforced canonical軌道の閉包 | 23 |
| \(C_T\) | constant-modular最高triple-lineのforced canonical軌道の閉包 | 22 |
| \(C_U\) | unimodular \(P_*\) のactual軌道の閉包 | 22 |

**全支持命題。** 採録済み全source定理の下で

\[
 \mathcal P_{\rm red}=K\cup C_G\cup C_N\cup C_F\cup C_C\cup C_D\cup C_T\cup C_U. \tag{1}
\]

従って全空間成分数を \(n\) とすると

\[
 142\le n\le149. \tag{2}
\]

台帳142個以外の成分は射影次元23以下である。七候補が相異なるとも、全空間成分とも、
既知成分から独立とも主張しない。従って149は等式ではない。
一般点の親、候補間の包含、全nilpotent厚み及び元の分類主定理は別問題として残る。
特に(1)は七候補の既知親への帰属又は成分性を要求しない。

## 2. actual因子の最上位場合分けは尽くされる

一般rank四の点では \(D=\operatorname{Pf}(\pi)\) は零でないactual五次因子である。
以下では被約化前のactual式を使い、その素因子と重複度で場合分けする。

まず重複超平面がある場合は§4へ送る。
重複超平面が無く、Dが非被約なら、重複する素因子の次数は二しかあり得ない。
次数三以上の二乗は五次を超えるためである。この場合は正確に
\(D=H Q^2\)、Q既約、H単純であり、§3.2が全配置を扱う。
残る被約Dの素因子次数分割は、次の七行で尽くされる。

| 被約actual Dの分割 | 全点の終端 | 使用する全source定理 |
|---|---|---|
| 5、正規 | K=T39 | `papers/normal_actual_divisor_recognition_ja.tex` の主定理。Dredだけの正規性へ拡張しない |
| 5、非正規 | K又は空 | `IRREDUCIBLE_NONNORMAL_QUINTIC_WHOLE_POINTWISE_COVERAGE.md`。導手曲面degree1..6、rank零/二、全plane位数の点別合成 |
| 1+4 | K又は \(C_G\) | §3.1の全cofactor/全marked-quartic源 |
| 1+1+3 | K又は \(C_N\cup C_F\cup C_C\) | §3.3の非plane十一型とplane十三配置 |
| 1+2+2、Q相異 | K又は空 | `TWO_QUADRICS_ALL_POINTWISE_KNOWN_COVERAGE.md` と独立全点合成査読及びscope addendum |
| 2+3 | K又は空 | §3.4のQ rank3/4/5全source |
| 1+1+1+2、及び1+1+1+1+1 | K又は空 | 三H＋Q dependent/independent全支持、及び五H全支持 |

これは素因子次数の表なので、rank二以下で分解するquadricを既約Qとして重複して残さない。
その場合は同表のより細かい素因子行、又は重複H行へ戻る。
以後のKへの矢印は全点包含である。旧G表の「像≤24」のセルをそのままKと読んでいない。

## 3. 単純因子を持つ全sourceの終端

### 3.1 H＋既約quarticの全残余は一つ

`HYPERPLANE_QUARTIC_COFACTOR_ENTRANCE.md` の24半単純代表と十nil affine型、
全90係数源、全Jacobi及びtranslation余核を使う。半単純全支持は
`HYPERPLANE_QUARTIC_SEMISIMPLE_SUPPORT.md` と `HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_SUPPORT.md`
の全分岐／SNC実曲線によりK又は空。

| nil型 | actual quartic既約の全source終端 |
|---|---|
| regular、c image内 | 新22成分とT39、全二源 |
| regular、長さ4top | 空 |
| nil31、c image内 | quarticに一次因子があり空 |
| nil31、短top | 全MPS実曲線でK |
| nil31、長top | 新22成分とT40、全二源 |
| nil22、c image内 | base括弧f非零は重複Hへ戻る又は旧ESSENTIAL5。f零は全marked plane-quartic源 |
| nil22、長さ2top | 七源全て既成分又はMPS |
| nil211、c image内 | 全源MINIMAL_M |
| nil211、長さ1top（Ga9） | 下記五源中の一源だけが \(C_G\)、他はK又は空 |
| nil211、長さ2top | 全三translation orbitから半単純実曲線、K |

型六のmarked入口では、mu6/r2十一源とmu5/r1五源は全frame／全六磁気／全境界までK。
mu6/r1も落とさず、A6/A4+A2 tangentはGa9 power1/3へ、E6 tangentは空又はGa9 linear-wへ戻す。
Ga9の全五源はdistinct quadratic、非零constant、及びc零のpower1/power3/linear-wである。
前二源は既23成分、power1/3全sourceは既cuspidal24成分へ入る。
`NIL211_SHORT_TRANSLATION_REVIEW.md` §§5–7はc零の各全11源の飽和閉包が
canonical23軌道閉包そのものに等しいことまで証明している。
従ってlinear-wの全11源を残しても一つの \(C_G\) だけで足りる。
全固定Dの二軌道と磁気unit/zeroは
`GA9_LINEAR_W_FIXED_DIVISOR_MAGNETIC_NORMAL_FORM.md` の全逆shearで覆われ、
磁気零軌道も同じunit軌道閉包に入る。source上界だけからこの包含を推測していない。

### 3.2 H＋同一Q二乗は全rank・全配置でK

`DOUBLE_QUADRIC_SUPPORT.md` はrankQ5と4の全配置を扱う。
rank5はH接触の既22軌道のみ、H非接触は空。
rank4はHがvertexを避ければ空、vertexを含む非接触は全beta零被覆、
接触は全一parameter族のT18と追加一軌道の実SNC pencilでK。
`RANK_THREE_QUADRIC_SUPPORT.md` は残rank3の全三配置を扱う。
横断は全beta零、接触primitiveは既22、余剰classの平方型は全source MPS、
混合型は既22、T(w)零の全11源は既24へ入る。
従ってQ既約の全rank3/4/5で、G-only又は未処理vertex終端は無い。

### 3.3 二H＋cubicのplane三源以外は全K

`TWO_HYPERPLANES_CUBIC_NONPLANE_WHOLE_COVERAGE.md` は全45係数と全低次項を残す。
非零半単純部は全coker・全共鳴を含めMPS。
nil部の十一affine型ではfactor又は不整合は空、planeへの全n6零補集合は十三配置へ戻り、
残非plane sourceは既22又は全actual SS/T06曲線でK。
特にtranslation rank二の全八源／七源、j21(h,x)六源、j21(0,h+x)六源とn6零補集合、
p0五源とp±1の両五源を戻す。source一般点の帰属だけを使わない。

planeの全十三配置は `TWO_LINES_PLANE_CUBIC_SUPPORT.md` の全Saito源で尽くされる。
九つの23次元成分とnode_two_node_tangents全六磁気のtoric collisionを除くと、
残るのはちょうど§1の三つの全A6磁気源である。
I非零のlambda正規形の親陰性をI零へ外挿しない。
三源の**元の全A6**を候補にするため、unit chartの外も(1)に含まれる。

### 3.4 Q＋cubic及び三H＋Q／五HはG上界ではなく全支持定理

Q rank5は `SMOOTH_QUADRIC_CUBIC_SUPPORT.md` の全MPS。
rank3/4は `DEGENERATE_QUADRIC_CUBIC_COFACTOR_REDUCTION.md` の全31代表を使う。
正則半単純三型、非正則半単純九型、rank4 nil五型、rank4 A零二型、
rank3 Anil3五型、rank3 A零七型は順に各REGULAR／SEMISIMPLE／Q4_NIL／
Q3_ANIL3／Q3_AZERO支持定理へ入る。
最後の六全磁気源も既23成分であり、31線形supersourceだけで止めていない。
従って全rank3/4/5のQ＋cubicはK又は空。

三相異H＋既約Qのdependent／independent全配置は両THREE_HYPERPLANES支持定理でK。
五相異Hは `FIVE_HYPERPLANE_SUPPORT_REVIEW.md` の線形式span5/4/3/2以下の全場合分けにより
T01又は空。この二行にも残るG-only sourceは無い。

## 4. 重複Hの全modular型：全sourceをK又は三軌道へ

**この節では標識 \(h^2\mid D\) を先に仮定する。** h-Casimir代表をh一上で
\(P=C+L+Q\) とすると、actual式は \(h^5\operatorname{Pf}P(x/h)\) だから
\(\operatorname{Pf}Q=0\)。modular場は \(Z=c+Ax\)、\(\operatorname{tr}A=0\) である。
全Jacobi、全C6/L24/Q60、全coker translationを保持する。
単純Hを持つ任意のrank一場へmultiple-Hの全17源定理を拡張しない。

| 全A／c型 | 全sourceの終端と厳密な補集合 |
|---|---|
| A可逆 | 全一次共鳴、対和零、退化C、非半単純を含む実MPS被覆でK |
| rankA3、半単純部非零 | simple-zero及び±a＋zero-Jordan2を含む全cokerでK |
| rankA3 nilJordan4 | c image外は全Pf零。image内は全五源が既五22成分 |
| rankA2半単純 | diag(lambda,−lambda,0,0)、c余核零／非零の全g定数・非定数・零と三共鳴がK |
| nil31、c image内 | n7単元は既DDH25、n7零はhalf／thirdの二つの全十源、既二25成分 |
| nil31、長さ3top | 全六源がK。r−1の両枝と全zero complementを含む |
| nil31、長さ1top | 元六source全てK。最後r2/3も全t単元・全d（d零含む）のT19曲線 |
| nil22、c image内 | 全11源がK。下表で全磁気と零補集合まで戻す |
| nil22、c image外 | 全八源がK。最後A1一軌道のT19 actual曲線も含む |
| rankA1、c image内 | 全CA、CB及びDDH25 |
| rankA1、c image外 | long全17源と全shortはK。各full25成分又は旧成分への全帰属、short→long actual曲線 |
| A零、c非零 | 全constant-modular sourceは \(K\cup C_D\cup C_T\)、§4.1 |
| A=c零 | 全unimodular sourceは \(K\cup C_U\)、§4.2 |

この表はtrace零四次行列の全Jordan型を尽くす。
rank二以下の非zero半単純型はdiag(lambda,−lambda,0,0)だけであり、
各nil型のtranslationは \(V/\operatorname{im}A\) 内の非零成分の最長block長で分類する。
異なる長さのtopやc=image外を勝手にc零へ置き換えていない。

nil31 image内のn7零二枝を旧 `MODULAR_NIL31_COMPONENT.md` §4だけで済ませない。
そこでは未被覆だったが、後続 `NIL31_IMAGE_ZERO_COMPONENTS.md`、
`NIL31_IMAGE_THIRD_COMPONENT.md` と独立 `NIL31_IMAGE_ZERO_COMPONENTS_REVIEW.md` が
half/third両全十源の成分性と全包含を証明している。
同様にnil22 image内は次の全11行で、旧「残stable」も後続全曲線により閉じる。

| 全Hessian源 | 既知終端 |
|---|---|
| stable_line_conic／stable_three_lines | 全六磁気T17／T03 |
| doubleA／doubleB／doubleC／doubleD | Ga9 cuspidal24／CPLUS／dependent三H＋Q24／DDH25 |
| tripleA／tripleB／tripleC | 既二nodal-cubic23／nil31-third25 |
| quadraticA／quadraticB | MPS40／nil31-half25 |

全190×90源からのHessian分類、GL2／C*全逆、任意F0六係数は
`NIL22_IMAGE_ELEVEN_WHOLE_POINTWISE_COVERAGE.md` の合成に保持される。
自由fiberを引くG上界をこの包含の代わりにしていない。
別横断条件[J,Q]=0の旧G-only sourceは、同じ全PをZ=c+Axで分類した上表の部分集合である。
その旧横断条件を追加の未分類terminalとして二重に残す必要は無い。

### 4.1 A零・c非零：全lower sourceは二forced orbitまで尽くされる

全constant場を \(-\partial_t\) へ規格化すれば、全sourceは厳密に

\[
 P=\partial_t\wedge Y+\operatorname{Jac}f,\qquad
 \operatorname{div}Y=1,\quad Yf=f,\quad \deg Y\le2,\quad\deg f\le3.
\]

fの加法定数もYf=fで固定する。degree二以下は全finite reduction、NEW25及び
sourceI/IIIのM11/B1 actual曲線でK。
f3被約は全nonconcurrent／concurrent実曲線でK。
f3非被約でもactual cubic自体が非被約なら全highest／lower sourceを
multiple-Hを再選択するrank一全17／short定理へ送りK。
残るf3非被約はGL3でx²y又はx³であり、二つの新全finite reductionが全lowerを扱う。

double-lineではf2 modxのrank／zero補集合、C単元／零、全f1/f0比較から、
affine、triangular、actual nonreduced、又はp*d0単元の一forced sourceへ帰着する。
後者の全四磁気parameterは全域のactual shearで三つを消し、残l単元を規格化する。
l零も同じ多項式全sourceの閉境界にあり、唯一canonical orbit閉包 \(C_D\) へ入る。

triple-lineではfull七syzygy quotientをQ3非cube／零／cube及びQ2 rank2/1/0で分ける。
前二群の全lower比較はaffine／triangular／actual cube又は空へ戻る。
Q3 cubeの全c単元／零、b単元／零、D単元／零を比較し、唯一bB単元の全四parameter sourceを残す。
その全actual f-automorphismとt-shearは分母を増やさず全parameterを同じcanonical点へ移す。
従って全sourceは \(K\cup C_T\)。
根拠はclaim1133–1135及び二つのFINITE_REDUCTION_SEPARATE_REVIEWであり、
全sourceの次元上界又は一代表のJacobi成功だけを結論にしていない。

### 4.2 A=c零：全h五乗と非飽和thetaを含め一軌道まで

全非退化closed二形式 \(\alpha=J+\beta+\gamma\) は

\[
 J^2\ne0,\quad d\beta=d\gamma=0,\quad J\wedge\beta=0,\quad
 2J\wedge\gamma+\beta^2=0,\quad\beta\wedge\gamma=\gamma^2=0
\]

を満たす。\(\theta=\iota_E\gamma/4\)、\(d\theta=\gamma\) を使う。
最高部の六foliation閉包だけから未知のC/L fiberを被覆したとはしない。
全低次延長に適用できる終端は次の通りである。

| 最高部の全閉支持 | 全C/L fiberの終端 |
|---|---|
| gamma零 | 全constant-Hamiltonian源CA |
| gamma factor付き | Grassmannian二rulingからlinear-Casimir又はfullpullback、双方K |
| R13全閉包 | proper linear-Casimir incidence、全affine Hamiltonian membership、K |
| R22全閉包 | 全実pencil。primitive regularはT16又は空、factorは上行、common-kernelはfullpullback、残Kroneckerは空 |
| S23全閉包 | proper constant-vector horizontal incidence、全pullback theorem、K |
| L1111/L112全閉包 | proper quartic integrating-factor incidence、全pole源はK又はP*一軌道 |
| E3全閉包 | nonzero nilpotent tangentのproper incidence。全四Jordanが空、linear-Casimir、closed L112又はpullbackへ戻る |

fullpullbackではw dwを省略せず、全C/L解を保つ。
Omega非零は全linear-Casimirへ、Omega零は非critical quadratic Hである。
rankH2が0/1なら既linear/constant-Hamiltonian、rank二ならH=z+xy、
全Bに対する \(B_e=B+eH\partial_z\) は全constant-modular sourceIIへ入るので閉性でK。
これは最高部のgeneric fibreだけを動かす議論ではない。

**非飽和thetaも別に尽くす。** 共通factorを持つthetaへ外部の飽和foliation六閉包定理を
直接使わない。theta=F_e*eta、e1/2/3を分ける。
この段落の四次式F4はthetaのintegrating factorであり、actual五次Dとは異なる。
e2ではeta=x dy−y dx型となりF4=xyF2のlog incidence。
e1ではd etaのlinear decomposable二形式の二rulingを使う。
共通一次方向ならeta=(x dG2−2G2dx)/3、F4=F1*x*G2のlog incidence。
三変数supportならsl3全linear像の稠密相異固有値でF4=F1*xyzを得て、
quartic integrating-factorのproper incidenceにより全Jordan／零residueも同じ最高部入口へ入れる。
これは**最高部incidence**の閉性であり、その後に全C/L log定理を適用する。
e3constant etaのEuler水平性からeta零、theta零もgamma零へ戻る。

log-pole最後のaligned-zero源の全残式とunit-Pfはn0=n3=a0零、n1*n2単元を強制する。
全translation／GL4 inverseを戻し、actual diagonalとtensor scaleでn1=n2一へ移す。
従って残る全sourceは一つの \(P_*\) 軌道。
これがclaim1129–1132及び `UNIMODULAR_H5_WHOLE_ONE_ORBIT_SEPARATE_REVIEW.md` の全scopeである。
元thetaに共通factorがあるときの穴、w dw補集合、E3の全境界を残していない。

## 5. 七候補の既約性と全像の次元上界

### 5.1 Ga9及びplane三源

Ga9の代表は `GA9_LINEAR_W_FIXED_DIVISOR_MAGNETIC_NORMAL_FORM.md` の

\[
 P_R=\partial_x\wedge U+\partial_z\wedge V+
       (-hx+2yz+R)\partial_x\wedge\partial_z,\quad R=4w^2,
\]

\(U=(4h^2-y^2)\partial_y+(16hy-4wy)\partial_w\)、
\(V=-hy\partial_y+(-4hw+4y^2)\partial_w\) とする。
固定D全R6のactual shearは五磁気係数を除き、残r w²のr単元はr四へ、r零は同閉境界へ入れる。
さらに \(y\partial_x+h\partial_z\) はactual Ga対称性なので
\(\dim C_G\le\dim\operatorname{PGL}_5-1=23\)。
ここでも元全11源の飽和が同じ閉包であることは§3.1の全source定理を使う。

plane三候補は、それぞれ固定した全磁気A6からPGL5×A6を作用させた像閉包である。
全磁気でPfaffianは同じ非零Dなので、これは既約domain上の全域の射影像であり既約。
次元を24+6と粗く数えず、実退化式写像の全fiberを使う。

固定した被約free plane quintic Dについて、Saito行列式のCramer割算は
AR0=AR1零、AR2が二次元であることを与える。
canonical liftのbase-base blockはprimitiveなdegree四のgradientを殺すdegree二blockなので零。
mixed blockのfiber次数一／二はAR1／AR0に入るので零。
従って固定Dの全rank四tensorは、二本の二次Saito場の定数GL2 frame四係数と
全磁気R6で尽くされる。Jacobiはfiber-fiberの線形partを一意に決める。
非零tensor scalarの一次元を除き、固定projective D fiberは高々9次元。
DのPGL5軌道は、三次元essential subspaceのGrassmannian6とPGL3軌道8から高々14次元。
従って各候補について

\[
 \dim C_N,\ \dim C_F,\ \dim C_C\le14+9=23. \tag{3}
\]

この全fixed-D source補題は既約quinticだけの分類を可約Dへ流用したものではない。
被約DのCramer/Saito自由基底及び元の十三配置定理を使う。
各配置が単一のprojective配置軌道であることもその有限配置定理に含まれる。
I零の全source、磁気零、frame entry零をこの議論から捨てない。

### 5.2 constant-modular二軌道

\(C_D\) のcanonical代表は座標(t,x,y,r)、pair順(tx,ty,tr,xy,xr,yr)で

\[
 (1+x-2xy,\ 4y^2-y,\ (1-6y)r+x,\ 0,\ -x^2,\ 1+2xy).
\]

\(C_T\) の代表は同じ座標で

\[
 P=\partial_t\wedge Y+\operatorname{Jac}(x^3+3xy-3y/2+r),
\]

\(Y=(x/3+y)\partial_x+(-xy+y/6-r/3)\partial_y+
(-3y^2+xr-5y/4+r/2)\partial_r\)。
両者とも全係数がtに依らず、actual translation \(\partial_t\) がtensorを固定する。
従って各軌道閉包は既約で、少なくとも次元上界24−1=23を持つ。
さらにtripleのHamiltonian xの場は

\[
 V=-(x/3+y)\partial_t+\partial_y+(-3x+3/2)\partial_r.
\]

これはaffine線形なのでactual PGL5内のPoisson対称性であり、
\(\partial_t\) と独立かつ可換である。実際のflowも
\((t,x,y,r)\mapsto(t-s(x/3+y)-s^2/2,x,y+s,r+s(-3x+3/2))\)
である。従って \(\dim C_T\le24-2=22\)。
七候補の上界は正確に五つが23以下、二つが22以下となる。
固定constant modular場の源次元やprimitive beta像の上界を使っていない。

### 5.3 unimodular残軌道には二つの独立actual Ga対称性がある

座標(x,y,z,w)、pair順(xy,xz,xw,yz,yw,zw)で

\[
 P_*=(x,-x^2,xz-1,xy-1,-yz+w-y,-xw+xy+y+z).
\]

二つのaffine場を

\[
 M=\partial_y-x\partial_z+z\partial_w,\qquad
 T=\partial_z+y\partial_w
\]

とする。\([M,T]=0\)。全六括弧に対するLie微分は、非自明な項を並べると

\[
\begin{array}{c|ccc}
 &xw&yz&yw\ \hline
 \mathcal L_M P_*&-x^2+x^2&x-x&(xy-1)-(xy-1)\\
 \mathcal L_T P_*&x-x&0&-y+y
\end{array}
\]

であり、zw成分はそれぞれ \((1-xz)+(xz-1)=0\)、
\((-xy+1)+(xy-1)=0\)。xy/xz成分も零である。
従って両場は厳密なPoisson対称性で、tensor scalarやEuler gaugeを追加する必要が無い。
その有限flowは

\[
 \exp(sM):(x,y,z,w)\mapsto(x,y+s,z-sx,w+sz-s^2x/2),
\]

\[
 \exp(tT):(x,y,z,w)\mapsto(x,y,z+t,w+ty)
\]

というactual可逆affine変換である。二場はPGL5内で独立、二flowは可換なので
stabilizerの次元は少なくとも二。
従って \(C_U\) は既約で

\[
 \dim C_U\le24-2=22. \tag{4}
\]

以上の上界は全PGL像の上界である。固定betaだけの像を全像と呼び換えていない。

## 6. rank二の全点と有限閉和の論理

§2–4によりrank四の全点が右辺Zへ入る。§5よりZは有限個の既約閉集合の和なので閉じている。
従ってrank四集合の全閉包もZに含まれる。

`RANK_TWO_SUPPORT.md` の正確な全点定理は

\[
 \mathcal P_{\le2}\subset K_{\mathrm{six\ rank2}}\cup
                         \overline{\mathcal P^{(4)}}.
\]

六成分はKの一部であり、第二項も上でZへ入れた。
同定の途中で使うtoric rank二及びAff32を六成分だけに無理に含める必要は無い。
特にrank二の既証明を「一般rank二成分が六個」だけへ弱めず、全点支持包含を使う。
これで(1)を得る。零tensorは射影空間の点ではない。affine coneを使う場合も原点は
既知閉coneに含まれ、追加の成分を作らない。

最後に \(Z_0\) を任意の全空間既約成分とする。
(1)の有限閉和に入る既約集合は、その一つの閉集合へ入る。
各右辺集合はPoisson空間の既約閉部分集合なので、\(Z_0\) の極大性から
その右辺集合と等しい。従って全成分は既知142又は七候補の一つであり、(2)が従う。
既知142以外なら(3)–(4)より次元23以下。
この推論は候補間の相違、局所smoothness又は全Jacobi idealの根基計算を必要としない。

## 7. 監査で特に補ったscopeと未解決の正確な意味

* normal D、全nonnormal integral D、H+Q1Q2、H+2Q、Q+C、三H+Q及び五Hは全点定理を使用した。
  どれもGだけのdimension terminalをKへ昇格していない。
* nil31 image内の旧n7零、nil22 image内のstable二源、shortnil r2/3、rank一全17/short、
  constant最高double/tripleの全lower、unimod nonsaturated thetaは明示的に後続定理へ戻した。
* simple-Hのrank一cofactorとmultiple-H全17定理を混同しない。
  plane三候補はそのsimple-H側に残し、全A6で定義した。
* 候補の閉包がrank二、別factor partition又は磁気零の境界を含んでも、右辺の重複であり漏れではない。
  七候補のうち幾つが既知成分に入るか、幾つが新成分かは未解決。
* 本命題は被約支持と成分数上界であり、schemeの被約性、generic transverse algebraの厚み、
  uniform denominator／splitting idealの新証明、exact149又は元の全分類完了を主張しない。

添付manifestは本証明、142/1135 accepted-input projection及び使用した全source定理の
不変proof headを固定する。各既存有限算術の再実行とは区別する。
古いnoteの「未被覆」は後続の全source証明を明示して更新して読むが、古いbytesは変更しない。
