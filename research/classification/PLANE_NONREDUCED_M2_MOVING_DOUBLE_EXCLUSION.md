# 非被約R・m二moving double rootは全P0空

2026-11-08、root純証明、独立読解前。新CAS・matrix・probe・jet無し。
S=P2、restriction rank二、f=a³又はa²b、actualD被約既約quintic、normal位数二、
非零double quadratic G2でprimitive normal rowのdegreeが一の場合を扱う。
degree零のnormaldouble、squarefree quadratic、m三は別残余。

## 全leading factorとseedのvertex

UFDと二gradingにより、全源で

    G2=gamma*h*ell², ell=A*u+B*v, gcd(A,B)=1,
    deg h=1, deg A=deg B=1, gamma≠0.

primitive degree一ならA,Bはindependent linear forms。
hのzeroはRに支持され、hはa又はb（f=a³ではaのみ）なのでq=f/hはdegree二。
P0は全normalweight零のprojective Poisson tensor、actualD(P0)=G2でrank四。
Q=Au+Bvはrank四quadric、vertex p=(A=B=0)∈S。
既受理DOUBLE_QUADRIC_SUPPORT §4のrank四vertex回避配置の全空性から **p∈h** が必要。
このseedへの必要条件をunknown higherの所属へ外挿していない。

## p∈qと、固定traceでの全connection一意性

qがhを含む場合はp∈hからp∈qが自動的に従う。
唯一自動でない場合はf=a²b,h=b,q=a²である。
もしp∈bかつa(p)≠0なら、p近傍のRはsmooth reduced line bのみで、betaはそこで単極。
tracefree部Gammaとroot characterは

    d ell+ell Gamma=(dlog q/2)*ell

を満たし、右辺はpでregularである。
local coordinates(b,t)でell=(b,t)とし、Gamma=(A db+C dt)/b+regularと書く。
b零でroot式を比べるとA,Cのresiduesはstrict upper triangular。
flatnessのb^(-2) db dt係数はC=±[A,C]なのでC mod b=0。
従ってdt componentにはpoleが無く、pでroot式を評価すると(0,1)=0、矛盾。
これはbに沿う局所単極だけを使い、別a²成分をreducedと仮定していない。
よって全sourceで **p∈h∩q** が必要。

tr beta=Tを固定する。二つのconnectionの差がactual primitive row ellを消しtrace零なら

    delta beta=omega*[[-AB,-B²],[A²,AB]].

matrix entriesのgcd一とf delta beta polynomialからf omegaはpolynomial one-form。
各coefficientのdegreeは零だからconstant。
Euler beta=IによりiE(f omega)=0、constant coefficient one-formの全coefficientsは零となる。
従ってdelta beta=0で、固定traceのconnectionは全Rで一意。

p∈qだからlinear auxiliary rowを補い、detM=qと出来る。
実際ellの二coefficientsが生成するidealにqが入るため、q=A*D−B*Cをlinear C,Dで表せる。
このとき

    beta=theta*I−M^(-1)dM, theta=(T+dlog q)/2

はroot式、trace T、flatness、Euler、f beta polynomialを全て満たす。
よって一意性からこれが元の全connectionである。
auxiliary rowにprimitivityは課さず、actual rootをmixするleftGL2も用いない。

## q=abの全source

この場合はf=a²b,h=a。p∈h=aなのでactual normal右GL2でell=a*u+t*v、t independent aへ戻せる。
全inverseは同じconstant GL2の逆。p∈qは既に満たす。
auxiliary completionを

    M=[[a,t],[0,b]]

と取る。beta21=0であり、正しいtracefree raw gaugeではnormal-v columnはv倍となる。
全normal K=k20u²+k11uv+k02v²は保持する。
actualG2のu²係数は正確に3f*k20だけ。
一方targetはgamma*h*a²なので、f=h*a*bで割れるためにはb|aが必要となる。
a,b independent、gamma≠0だから不可能。
tがbと比例する境界も、そのほかの全p位置も同じ非整除に含む。

## q=a²の全source

f=a³,h=a、又はf=a²b,h=bの全場合である。
p∈aからell=a*u+t*v、t independent aへactual normal GL2で戻す。
auxiliary completionは

    M=[[a,t],[0,a]],
    beta=(theta−dlog a)*I−[[0,d(t/a)],[0,0]].

従って全raw sourceは

    P_base,u=uX+vZ, P_base,v=vX,
    P_uv=k20u²+k11uv+k02v²

という同一diagonalの形となる。thetaの全irregular exact部分を除いていない。
actualG2のuv係数ではXとXのcross determinant及びEuler差が消え、Zはv²項にしか寄与しない。
よってこの係数は **3f*k11=3h*a²*k11**。
しかしtarget gamma*h*(a*u+t*v)²のuv係数は2gamma*h*a*t。
h*aをcancelして3a*k11=2gamma*tとなり、t independent aに反する。

これで全f=a³/a²b、全h、全moving degree一double normal rootのP0は空。
全fullsolutionはそのP0を持つ必要があるため、higher係数を一つもprobeせず全higher源も空。
constant root degree零、normal quadratic squarefree及び他multiplicityは結論しない。
