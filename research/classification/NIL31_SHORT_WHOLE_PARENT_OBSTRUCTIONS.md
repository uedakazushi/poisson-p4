# nil31短鎖 r=2/3：全 stabilizer・rankzero plane・marked first jet による親閉包の排除

2026-11-08、PRIMARY E、`/root/quartic_review` の純証明。
新 CAS、matrix、minor、profile、jet 又は数値点は実行しない。
target は `MODULAR_NIL31_SHORT_SUPPORT.md` §7 の canonical(15)、**td!=0**。
同節の全 affine inverse と二parameter torusによりこれは全 source の稠密22次元軌道を代表する。
全 E所属は未解決であり、以下の有限な全親排除だけを主張する。

## 1. stabilizer 全 Lie algebra と全 constant-conjugate 1PS の障害

既受理の orbit lower bound は `CLASSIFICATION_NIL31_SHORT_COMPONENTS_20261107.json`
の `r_two_thirds`、orbit/vertical/scale rank47、vertical24から projective orbit≥22。
元 integer source の47minorは保存素数で775950!=0、従って標数零でも非零。
全 canonical sourceとの実逆、同じ witness t1,d−58への移送は
`CLASSIFICATION_NIL31_SHORT_TWO_THIRDS_CANONICAL_20261107.json` に固定。
さらに `CLASSIFICATION_NIL31_SHORT_TWO_THIRDS_QUADRATIC_20261107.json`
は同じ点の QQ orbit/vertical/scale rank47、全 pivot rows と実行sourceを保存する。
新点や小matrixでこの orbit boundを再推測しない。

一方、全 target の二つの projective automorphisms は

    N=Y∂X+Z∂Y+h∂W, A=h∂X。

Nは affine Poisson の BV identity、Aは h-Casimir lift のX非依存から tensorを保つ。
[N,A]=0、PGL Lie内で独立。従って orbit≤24−2=22。
下界と合わせ orbit=22、stabilizer Lie は正確に C N+C A。

座標順序(X,Y,Z,W,h)で行をvectorfield成分、列をcoordinateとすると

    N=[[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],
       [0,0,0,0,1],[0,0,0,0,0]],
    A=[[0,0,0,0,1],[0,0,0,0,0],[0,0,0,0,0],
       [0,0,0,0,0],[0,0,0,0,0]]。

全 linear combination が strictly upper triangular、従って nilpotent。
この Lieには非零 semisimple PGL元が無い。
任意非trivial cocharacter λ:C*→PGL5 と固定 g,p に対し、
非零 projective limit lim λ(s)g·p は λ(C*)に固定される。
その Lie は非零semisimpleなので、target orbitへは入らない。
したがって **任意 source の fixed-conjugation後の純1PS limitでこのtdunit軌道へ達することは無い**。
moving preconjugation g(s)まで排除したものではない。

T37 zero Veronese+secantの提案familyも、合成後は
diag(s8,s4,s6,s8,s7)·Ctheta·T^-1 という固定conjugation後の純1PSである。
よって actualtensorpushを実行せず、この正のtarget目標を棄却できる。
zero-supportの幾何一致だけでtensor所属を推論しない。

## 2. target の projective rankzero reduced support

affine h1では PZW=2Z/3、PYW=(Y+ZW)/3、PXW=d+YW/3。
rankzeroなら Z=Y=0 とd=0を同時要求するので、tdunitでは外に点が無い。
h0では projectivezeroを rawP=0ではなく **E∧P=0**で判定する。
残る四cubicは

    Z{(Y²−2XZ)/3+2tZW},
    W{XZ/3−2Y²/3+2tZW}, −YZW, −Z²W。

W0では Z(Y²−2XZ)=0、Z0/W!=0ではY0。
従って reducedsupport は正確に

    C={h=W0,Y²−2XZ0},
    L1={h=W=Z0}, L0={h=Y=Z0}。

Cは滑らかなconic、L1はCにp=[X]=[1:0:0:0]で接し、L0はそのplane外のlineで同じpに会う。
特に rankzero P2は存在しない。
非被約zero scheme、各component multiplicity又はflatzero-familyのdegreeは本証明では未主張。

## 3. Ga9-distinct／nil22-long case5 の whole23親排除

いずれも既受理のdense23 actual orbitの閉包が全成分である。
accepted primitive tensorsそのものから、rankzero P2を次で確認出来る。

Ga9-distinctは coordinates(a,b,c,d,e)で plane b=d0。
normal rowsは零、残る括弧は

    Pac=48a²,Pae=32ac,Pce=−48ae+32c²。

aPce−cPae+ePac=0なので、そのplane全体でprojective tensorは零。
nil22-long case5は coordinates(h,x,y,z,w)で plane h=y0。
normal rowsは零、残りは

    Pxz=−2x²,Pxw=−243xz,Pzw=2wx−243z²。

xPzw−zPxw+wPxz=0。
rawtensorが零でなくても projective tensorが零であることを正しく用いた。

I={(Pi,S):S∈Gr(3,5), Pi|S=0} は閉incidence、Grassmannianはproper。
従って rankzero planeが存在するprojective Poisson tensorsの集合は閉。
両23 denseorbit、従って全component closureにこの条件が成り立つ。
targetには§2よりrankzero P2が無いので、**両whole23component closureから排除される**。
singular planeの本数をclosed条件として使ったのではない。

## 4. 二つの新 triple weak-wall23成分の全 marked plane の特殊化

対象はτ=4/21,1/12の全component Cτ。
`PLANE_TRIPLE_WEAK_WALL_COMPONENTS_SEPARATE_REVIEW.md` により、各々は以下のdense actual23orbitの閉包。
coordinates(a,b,c,z,w)、Pi=Qτ+∂c∧Wで

    Qab=−(3τ+2)az+w²,Qaw=−3z²,
    Qbz=(6τ+1)z²,Qbw=9τzw,Qaz=Qzw=0。

**完全な affine particular**（W0だけにしない）は

    τ4/21: W=(a²+7bw,−2ab,4az/7+w²/2,aw/7−3bz),
    τ1/12: W=(a²−20bw,−2ab,az/4−w²,−aw/2−3bz)。

唯一のsingular planeは S={z=w0}、actualD は globally I_S³に入る。
SはPoisson invariant、restriction zero cubicは −3a²b。
double line La={a0}とsimple line Lb={b0}をS内で標識する。

target∈Cτ と仮定し、generic pointがこのdense orbitにあるDVR curve Pi_sを取る。
actual S_s、La_s、Lb_s はproper flagvarietyで finitebasechange後延長出来る。
D_s∈I_(S_s)³、S_sのPoisson不変性もclosed equationsなのでspecialfiberへ残る。
特異schemeのflatnessは仮定しない。

target D=h²Z(hY−ZW) に mult≥3のplaneは次の二つだけ。

    S0={h=Z0}, S1={h=W0}。

実際、h外ではZとQのsimple成分の交差でもgenericmult≤2。
Qのsingular vertexは点なのでplane全体でQをdoubleに出来ない。
h内でZ,Wの両方がgenericunitならmult2。
mult≥3にはplane⊂Z0又はW0が必要で、上の二planeとなる。

そのrestriction zero cubicsはそれぞれ

    R_S0=−2Y²W/3, R_S1=Z(Y²−2XZ)/3。

両restriction tensorは非零。
relative planeをregularにtrivializeするとzero cubicの係数はregularでspecialfiberも非零。
doubleline*lineのcubic集合は(P2dual)²のproper multiplication imageなので閉。
R_S1はirreducibleconicを持ちこの像外、従ってlimit S=S0。
R_S0のdouble/simple factorsは相異なり一意なので

    La_0=L0(Y0), Lb_0=L1(W0)。                       (1)

marksがcoalesceして別のlineへ逃げる補集合も、このnonzero restriction cubicにより除かれた。

## 5. generic simple line の full firstjet matrices

parent Lbのa!=0 chart a1で、点p=(b,z,w)0,c=c0を取る。
元 complete particularから全 bracketsを一次まで展開すると

    {b,c}=3b−(3τ+2)c0 z,
    {c,z}=(3τ−1)z,{c,w}=(6τ−2)w、他の一次項は零。

省略してはいけないbw/w²等の particular termsも戻したが、これらの一次項には寄与しない。
isotropy Lieのderivedidealは D=span(db,dz,dw)、abelian dimension3。
ad_dc|D の matrix（columnsが像）は

    [[−3,0,0],[(3τ+2)c0,3τ−1,0],[0,0,6τ−2]]。

従って eigenvalues は λ(k,1,2)、λ=3τ−1!=0、
k=7（τ4/21）又は4（τ1/12）。c0によるoffdiagonalはtrace invariantsを変えない。

target generic L1のY1 chart、点(h,Z,W)0,X=x0では

    {h,X}=−h/3,{X,Z}=−Z/3,{X,W}=(2W−x0h)/3、
    他の一次項は零。

derivedideal D0=span(dh,dZ,dW)はabelian dimension3、ad_dX matrixは

    [[1/3,0,−x0/3],[0,−1/3,0],[0,0,2/3]]。

eigenvaluesは(1,−1,2)/3、trace非零、tr(A³)/(trA)³=1。
rootも両親・targetの全一次項と固有値を独立に手で再導出した。

## 6. derivedideal bundle と homogeneous trace identity の極限

(1)により、generic target p0∈L1\{p}をLb_s上のregular section p_sへlift出来る。
このsectionはgeneric La_s∩Lb_sを避けるので、parentの全generic firstjet formulaが適用される。
P_s(p_s)=0、従って cotangent rank4 bundle上にregularなisotropy Lie bracket族を得る。
source座標のc0がpoleを持つ場合も、target近傍のregularframeでのbracketはregularであり、
generic isomorphism classだけを上の計算から用いる。

derivedidealは bracket map Λ²L→L のimage。
genericとspecialでともにdimension3なのでunit3minorを持つconstant-rank mapとなり、
image Dはdirect-summand free rank3 bundle、L/Dはfree rank1。
generic Dのabelian性はpolynomial bracketidentityなのでwholebundleへ延長する。
quotientgeneratorを持ち上げるregular section uを選ぶ。
A=ad_u|D はregular3×3 matrix。
genericでは任意generator/座標/Poisson scalarについてAのeigenvaluesは γ(k,1,2)であるから

    (k+3)³ tr(A³)−(k³+9)(tr A)³=0。                 (2)

これは **generatorにhomogeneous degree三** のpolynomial identity。
uのgeneric対応scalarγがpoleを持つ場合も等式はK上で成り立ち、Aのregular性によりRへ延長する。
derivedideal内のlift変更はDがabelianなのでad_u|Dを変えない。
traceで割る、genericbasisをspecialfiberへ固定する、又はflatnessを仮定する操作は無い。

targetではu0はnonzero quotientgeneratorだからA0は上のtarget matrixのnonzero scalar倍。
従ってtr(A0)!=0、tr(A0³)=(tr A0)³。
(2)の左辺はtargetで

    { (k+3)³−(k³+9) }(tr A0)³

となる。k7なら係数648、k4なら270、いずれも非零。
これは矛盾。
よって **shortnil tdunit軌道は両weak-wall whole23component closureに入らない**。
一つの一般因子型やselectedjetが違うだけの比較ではなく、全proper plane/line flagと
constant-rank derivedidealの特殊化からwhole閉包を排除した。

## 7. 正確に残る E課題

本稿は pure1PSの全障害と、四つの指定23全成分からの非所属を与える。
marked multiple-H SS2の別全closure排除、固定SS12companionの別diagnosticsも不変に保持する。
全MPS23/24/25親、moving singularflags、全r2/3 sourceのE所属又は新成分性は未解決。
これらの負の全親証明を、任意knowncomponentからの非所属や分類完成と呼ばない。
