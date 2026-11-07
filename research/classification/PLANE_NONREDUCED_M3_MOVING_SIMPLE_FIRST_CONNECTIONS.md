# 非被約R・m三generic21・constant double root/moving simple rootの全finite入口

2026-11-07、producer純幾何。primitive degree(0,1)/(0,2)の全first connectionsを扱う。
新globalK/operator/G上界はまだ計算していない。constant(0,0)、moving-double(1,0)、
generic triple、m二、reduced Rへ結論を拡張しない。

## 1. 共通式

constant double rootをz、primitive moving simple rowをAz+Cwとする。
全leadingは `G3=h z²(Az+Cw)`、d=degA=degC=1又は2、degh=2−d。
gcd(A,C)=1で、h及びCのdivisor支持はRに含まれる。
lower triangular β=[[p,0],[r,q]] の全parallel/flat条件は

`p=dlog(f/(hC)), r=−d(A/C)+(A/C)(q−p), dq=0`。

qの全nonlog部分と、global fβのpolynomial条件を使う。
flatnessからglobalK存在を推論せず、後続finite batchで全三K/P0Jacobiを保つ。

## 2. degree(0,1)：全四whole line＋十五fixed

R2Lではf=a²b。h=a又はb、C=a又はb。
AはCのline上でR他lineとの交点又はそれ以外という二orbitなので、
C=aにはA=b,c、C=bにはA=a,cで尽くされる。
全closed qは

`q=(dlogf)/3+αdlog(a/b)+γd(b/a)+δd(c/a)`。

C=aではf*rのa最高poleからγ=δ零、α全C。
h二選択×A二選択で四whole lineを保持する。
C=bではb最高poleから、h=aならα=4/3、h=bならα=7/3。
全γ,δはglobalに許される。
actual root-row stabilizerのc shear/diagonals/normal Borelで
δ1γ0、δ0γ1、両零の全三strataへ入る。
したがってh二×A二×三で十二fixed。
逆のnormal/base変換は `PLANE_NONREDUCED_M3_MOVING_DOUBLE_FIRST_CONNECTIONS.md` §5と同じで、
p/log residuesを保ちqのexact coefficientsだけを戻す。

R3Lではf=a³、h=C=a、A=b。p=dloga。
全q=dloga+d(P2/a²)のうちglobal f*r条件はa|P2を強制し、
q=dloga+γd(b/a)+δd(c/a)。同じ三strataで三fixed。
以上四whole line＋十五fixedで全degree(0,1)源を覆う。

## 3. degree(0,2)：R2Lの全closed-source制限

h定数、f=a²b。Cはa²、ab又はb²。
Aは任意quadraticでgcd(A,C)=1。qは上記全三parameter closed connection。

C=a²ではp=dlogb。
f*rのa最高poleからγ=δ零、次のpoleからα=−8/3。
従ってq=−2dloga+3dlogbで、Aはa非可除の任意quadratic。

C=abではp=dloga。
b poleからα=4/3、a poleからγ=δ零。
従ってq=2dloga−dlogbで、Aはa,bのいずれでも可除でない任意quadratic。

C=b²ではp=2dloga−dlogb。
b最高poleからα=10/3。次の全coefficientをb零へ制限すると

`a∂c(A|b0)=δ(A|b0), a∂a(A|b0)=2(A|b0)`

となる。非零polynomial A|b0について前式はδ零を強制し、後式と合わせA|b0=k a²、k非零。
さらにdb coefficientは∂bA|b0=kγaなので、
A=k(a²+γab)+m b²。
normal rowのrelative scaleでk一、normal Borelでm零へ戻す。
q=4dloga−3dlogb+γd(b/a)を保持し、actual a/b ratioでγ零/一へ。
この二fixedが全C=b²である。

これらは同じ閉formをlogと仮定したものではなく、元全nonlog directionsをglobal pole条件で排除又は保持した。

## 4. primitive quadratic Aの全normalizer orbits

normal Borel w→w+κzはAのC倍を変えるため、Aをmodulo Cで考えられる。
base a/b diagonalsとc→sc+ua+vbはRを保つ。

C=abでは、Aのc²係数非零ならcomplete square後
A=c²+εa²+ηb²、ε,η∈{0,1}の四型。
c²零/c-linear両a,b係数非零ならc shearsでa²,b²を消しA=c(a+b)。
c-linearがaだけならprimitive性によりA=ac+b²、bだけならA=bc+a²。
c-linearも零ならprimitive性でa²/b²係数双方非零、A=a²+b²。
合計八fixed。

C=a²ではc²係数非零ならA=c²+εab+ηb²の四型。
c²零/c-linearのb係数非零ならab/b²を消せ、A=bc又はc(a+b)の二型。
b係数零/a係数非零ならprimitive性よりA=ac+b²の一型。
c-linear零ならprimitive性からb²係数非零、A=b²又はab+b²の二型。
合計九fixed。

全caseはc²、c-linear各coefficientの零/非零で尽くし、
normal/base非零scalesとminus shearsでactual inverseがある。
primitive性に反するA可除a/b枝はdegreeを落とす元UFD factorizationの補集合であり、本flagへ含めない。
R2Lの全degree(0,2)は九＋八＋二＝十九fixed。

## 5. degree(0,2)：R3L

f=a³、C=a²、p=dloga。
q=dloga+d(P2/a²)を全て戻すと、global f*rはa|P2、続いてP2∈C*a²を強制する。
従ってq=pであり、r=−d(A/a²)。
primitive Aはaffine quadratic A/a²のquadratic part非零と同値。
affine GL2とtranslations、modulo constantにより
A=bc、b²、b²+acの三型。
rank二quadratic、rank一でkernel-linear零/非零の全caseである。
R2L十九と合わせdegree(0,2)は全二十二fixedとなる。

## 6. 必要な次の検証

各modelはrawω=3fβ−dfからtracefree linear mixed sourceへ戻し、
全K三未知のactualD2零、P0全350Jacobi、fullP1三十/P2九、
actualD3が全line C*h z²(Az+Cw)へ入る必要条件を確認する。
finite源の列挙だけでglobal Poisson lift又はG閉鎖を宣言しない。
全leading amplitudeを未知のまま保ち、zero leadingはm四以上の別義務へ送る。

安全なfixedplane first-symbol dimension upperは、degree(0,1)でmoving-double同配置と同じ：
C=a whole linesはA=b/cで9/10、固定α壁は一減。
C=b exact三strataはA=aで10/9/8、A=cで11/10/9、R3Lで9/8/7。
degree(0,2)のR2L C=a²/abは全quadratic root-row sourceを保持したB=12が安全。
C=b²のγ零/非零はB=8/9。
R3L A=bc,b²,b²+acはそれぞれB=10/8/9。
これらはmarked parameter sourceの上界で、canonical点のorbit dimensionへ特殊化した値ではない。

現時点の成果は全finite firstconnection分類のみで、全G/Eの完了とは異なる。
