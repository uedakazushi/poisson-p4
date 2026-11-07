# 非被約zero cubic・m二split normal roots：全極予算と三flagsの排除

2026-11-07、主義務G、担当 `/root/quartic_review`。**今回の純数学的証明**。
rank二Poisson plane S=P2、非零zero cubic f、actualD被約既約quintic、normal位数二、
generic squarefree normal quadratic G2を仮定する。
本稿ではf=a³又はa²b（a,bは相異なるlinear forms）を対象とする。
新しいmatrix/Jacobi probe・minor・GB・jetは実行しない。
全G/E、normal double quadratic又は残flagsのfullhigherを完了したとは主張しない。

**命題。** root degree flags(0,3),(1,2)の全connectionは空。
flag(1,1)はflat connection自体が存在しても、全global P0 consistencyで空。
残る(0,0),(0,1),(0,2)の全connectionは下表のirregular exact formsを含む源に尽くされる。

## 1. 全nonreduced入口とpole budget

既 `PLANE_M2_SPLIT_ROOT_MATRIX_ENTRANCE.md` のHamiltonian associated-graded証明は、
Rの被約性を使わず `Delta(G2)=c f²`、c非零を与える。
UFD及び二gradingにより

`G2=h ell1 ell2, (ell1,ell2)^t=M(u,v)^t, h detM=f`

と出来る（非零constantsはhへ吸収）。各root rowはprimitiveであり、共通零点を許す。
row degrees d1≤d2とdeg h=3−d1−d2から、全flagsは
00,01,02,03,11,12である。
normal double-leading Delta零はこの証明の別補集合として残す。

全conormal connection betaは、f betaの各one-form係数がdegree二のpolynomial、beta(E)=I、
flatness d beta−beta wedge beta零を満たす。
prime pのfにおけるmultiplicity nに対しbetaのpole orderは高々n。
primitive rowのunit成分からtheta_iのpole orderも高々nと**先に**証明し、

`dM+M beta=diag(theta1,theta2)M`, `d theta_i=0`, `iE theta_i=di+1`

を得る。
M inverseだけからthetaのpole orderを推測しない。

eta=theta1−theta2、m=ord_p(detM)とする。m非零ならM modpはrank一、両rowとも非零。

`C=adj(M)diag(1,0)M`

はmodpで非零rank一、trace零かつnilpotentであり、少なくとも一つのoff-diagonal成分がunit。
これは両rowのprimitivityを使う。
`beta=theta2 I+eta*C/detM−M^-1dM` のそのoff-diagonal成分を比較する。
etaにorder k>0のpoleがあればmiddle termはorder m+k。
M^-1dMは高々m、scalar theta2はoff-diagonalに寄与せず、betaは高々n。
従って **k≤n−m**。
m零でもtheta_iの元boundからk≤n。
全divisor等式h detM=fにより、全primeで

**`pole_divisor(eta) <= div(h)`**

を得る。pがhとdetMを共に割る、又はm二／三のcollisionもこの式に含む。
被約性又はcoprimalityを使っていない。

## 2. flags03/12はnonreducedでも全空

h constantではetaの全divisorial polesが消える。
原点及びroot basepointのcodimension二に新しいrational coefficient poleは無い。
etaはEuler weight零、すなわちcoefficient degree−1なので、polynomial one-formなら零。
一方iE eta=d1−d2であり、03/12では非零。
従って両wholeflagsの全connection、従って全Poisson sourceが空。

hがlinearならetaはそのline沿い高々simple。
closednessからeta=c dlog h、Eulerからc=d1−d2となる。
特にflag11では **eta零**、theta1=theta2=theta。
この結論はfがa³/a²bでも成立する。

## 3. flag11・detがsplit distinct二線の全枝

eta零なので全betaは `theta I−M^-1dM`。
scalar thetaのirregular exact部分を捨てない。
f非被約なのでdetMはab又はa²だけ。
abの場合は既split-matrix証明のactual正常形がそのまま使える。

image dimension三ではannihilator rank一。
primitive両row条件からannihilatorのleft成分は共に非零なので、
合法なrow rescalingとnormal右GL2だけで

`M=[[x,y],[x,z]], detM=x(z−y)`

へ戻る。rootをmixするleft GL2を捨てていない。
M^-1dMはupper triangular、従ってbeta21零。
canonical tracefree mixed gaugeではnormal-v columnはv倍なので、
actualG2のu²係数は必ず **3f k20**。
しかしtargetの同係数はh x²であり、f=h x(z−y)で割れない。
hがx又はz−yと重なる場合でもhをcancelして同じ矛盾になる。

image dimension二では全regular pencilを
`M=U diag(x,y)V` と書く。Vだけをactual normal右変換で除き、**全left Uを保持する**。
両row primitiveからUの全entriesが非零。
M^-1dMはdiag(dlogx,dlogy)で、theta scalarによりUは接続式からだけ消える。
actualu²係数3f k20に対しtargetはh U11U21 x²、f=h detU xyなので再び非整除。
image一はprimitive linear rowsに反する。
従ってsplit determinantのflag11は、hとdetの重複を含め全P0が空である。

## 4. flag11・double-line determinantの全left-U族

detM=a²とする。Mのlinear image W⊂Mat2のdimensionは高々二。
実際detはMat2上のnondegenerate quadratic form。W三次元ならradicalはdimension高々一なので、
restrictionのrankは二以上となり、rank一のa²に反する。
W一ならM=a Cで両rowがdegree一共通因子を持つのでprimitiveでない。

従ってW二次元のregular pencilであり、invertible memberを一に合わせた残matrixは
nonzero nilpotent Jordan blockとなる。
base coordinates a,t、constant invertible matrices U,Vにより全sourceを

`M=U [[a,t],[0,a]] V`

と書ける。tはaとindependentなlinear form。
**Uをroot同値として除去しない。** U=[[p,q],[r,s]]、detU非零を全て保持する。
normal右GL2でVだけを除くと、root rowsは

`ell1=pa*u+(pt+qa)*v`, `ell2=ra*u+(rt+sa)*v`。

primitivityはp,r非零を強制する。
p又はr零なら対応rowはaの共通因子を持つからである。

eta零によりleft Uはbetaからのみ消え、

`beta=(theta−dlog a)I−[[0,d(t/a)],[0,0]]`。

特に両diagonal entriesは等しく、beta21零。
canonical tracefree mixed gaugeでは

`P_base,u=uX+vZ`, `P_base,v=vX`

となる。全normal bracket `K=k20u²+k11uv+k02v²` は保つ。
Euler-signed Pfaffianのuv係数は **3f k11**。
det(E,X,X)は零、normal Euler項のuv寄与も同一Xの差として相殺し、Zはv²にしか寄与しない。
従ってこの係数はfで割れる。

他方target h ell1 ell2のuv係数は

`h*a*(2pr*t+(ps+qr)*a)`。

f=h detU a²で割れるにはaが2pr*tを割る必要があるが、tはaとindependent、pr非零なので不可能。
これでdouble-line determinantの**全left-U族も全P0が空**。
scalar thetaのirregular性、h=a又は他line、root basepointの位置をgenericにしていない。
新しいP1/P2又はchosen minorを計算する前にglobal consistencyだけで閉じた。

## 5. 残flagsの全connection表：exact部分を残す

以下の記述は接続源の必要十分な表示であって、P0又はfullhigher sourceの存在定理ではない。
任意のclosed rational form theta with pole divisor≤f、Euler value cの全体をK_f(c)とする。
independent linear coordinatesを(a,b,c0)とすれば全formsは

* f=a³：`c dlog a+d(P2/a²)`、P2は全quadratic modulo C a²（五exact係数）。
* f=a²b：`lambda dlog a+(c−lambda)dlog b+d(L/a)`、Lは全linear modulo C a（二exact係数）。

例えばa一chartでlog residueを引いたclosed formをpolynomial積分し、
a沿いのpole order上界からpotentialのdegreeを二又は一へ抑えれば完全性を得る。
全Euler-degree条件を保つ。nonlog exact部分をzeroへ固定しない。

flag00はconstant rootsu,v、h=f、beta=diag(theta1,theta2)、theta_i∈K_f(1)が全源。

flags0d、d一又は二では

`M=[[1,0],[A_d,B_d]], hB_d=f, gcd(A_d,B_d)=1`。

A_dは全homogeneous polynomial、B_dは下表の全component product。
theta1は任意K_f(1)、theta2=theta1−etaとすると、全betaは

`beta11=theta1`, `beta12=0`,
`beta22=theta2−dlog B_d`,
`beta21=−(dA_d+A_d eta)/B_d`。

これは全closed diagonal connectionをactual rational gaugeで戻した式でありflat。
Euler beta=I、f betaの全coefficientsはdegree二polynomial。
特にoff-diagonalのdenominator B_dはf=hB_dで正当に消え、h etaはregularなので別poleを増やさない。

| f | flag | B_d | h | 全eta |
|---|---|---|---|---|
|a³|01|a|a²|−dlog a+d(L/a)|
|a³|02|a²|a|−2dlog a|
|a²b|01|a|ab|r dlog a+s dlog b、r+s=−1|
|a²b|01|b|a²|−dlog a+d(L/a)|
|a²b|02|a²|b|−2dlog b|
|a²b|02|ab|a|−2dlog a|

Lは全linear modulo C a。
A_d/B_dの共通因子境界はprimitive factorizationでhへ戻した別flagであり、無断に捨てない。
base/normal実変換のinverse、二rootの交換、nonzero amplitudeもこの表の全sourceへ戻せる。

## 6. 次の必要義務と保存範囲

残00/01/02のP0 global consistencyは、canonical tracefree mixed field gaugeで全三係数Kを残し、
全Euler-Pfaffian二次normal係数がh ell1 ell2に一致し、元projective Jacobiが零であることを要求する。
flat connectionだけから存在を推論しない。
fullhigherは依然全P1三十・P2九と全quadratic Jacobiを保持する。
本稿はこの有限化されたirregular connection表を次の構造数学の入口として保存する。

今回証明したのはnonreduced f・m二SFの03/12全connection空、11全P0空、
及び残00/01/02の完全connection poolである。
残flagsのG上界/E帰属、Delta零のdouble normal leading、nonreducedm三、全G/Eは未証明。
独立別作者reviewは今後のreceiptを要する。
