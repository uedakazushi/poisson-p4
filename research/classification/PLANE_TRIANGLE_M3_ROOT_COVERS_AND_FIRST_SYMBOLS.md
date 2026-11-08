# Triangle零曲線・plane位数三：全root coversを五diagonal first symbolsへ

2026-11-08、G、producer plane_collision_review。純幾何証明、独立review前。
新source probe、minor、GB、jet又はCASを実行していない。
S=P2、R=(xyz=0)は三つの非concurrent lines、πS rank二、actualD被約既約quintic、normal位数三を仮定する。
結論は全first-symbolの有限必要還元、generic21と下記SF枝の全空性である。
最後の五full higher源のG上界はまだ証明していない。

## 1. 全root typesとKummer characters

U=(C*)²なのでSF root monodromyは1,C2,C3だけ。
C2のnontrivial characterは、三line meridiansのrelationから二つのcomponentsだけでbranchする。
base permutationでbranch pairをx,yと取る。
C3のnontrivial characterは、inverseとbase permutationsを除き
weights(1,1,1)、又は(1,2,0)だけ。
前者のcompactificationはw³=xyzというnormal cubic 3A2 surface、
後者はnormalization P(1,1,3)（有限map[s:t:z]↦[s³:t³:z]）である。
scalar GL2 monodromyをroot permutationと混同せず、全root coversをこの有限表へ含める。

generic21は全degree flags00/01/02/10、tripleはconstant root、split SFは000/001/002/011。
その全primitive degree証明とzero complementsは `PLANE_REDUCED_CUBIC_M3_ROOT_GEOMETRY.md` に固定する。
本稿でnormal反応を検査する前にchosen root profileを仮定していない。

## 2. Toric characterと全actual coefficient・curl式

ordinary diagonal charactersを
`p=sum p_i dlog x_i, q=sum q_i dlog x_i`, sum p_i=sum q_i=1とする。
f=xyz、P_fのbase bracketsは{ x,y }=xy、{ x,z }=−xz、{ y,z }=yz。
tracefree mixed fieldsX,Yのeigenvaluesは

```
X=(p_y−p_z, p_z−p_x, p_x−p_y),
Y=(q_y−q_z, q_z−q_x, q_x−q_y).
```

normal K11=κはD2零から

`κ=(1,1,1) dot(p cross q)`。

この符号はD2のuv coefficient `3fκ−det(E,X,Y)`を比較したもの。
off-diagonal fieldがあってもこのuv式は変わらない。
constant root uの全bad vector Cは

`q_i C_i=0 (all i), (2κI+4DY−DX)C=0, C≠0`。          (1)

constant root vなら対称に
`p_i C'_i=0, (-2κI+4DX−DY)C'=0,C'≠0`。
前半はfullD3係数、後半は全weight一Jacobiのcurl。
off-diagonalと他の全P1/P2はこの二identityへ寄与しない。
rootでないparallel axisでもD3のv³ coefficientは3f q(C)、u³ coefficientはsignを除き3f p(C')。
従ってdiagonal sourceの両axisのpure cubic coefficientsは常にspan(xy,xz,yz)へ入る。
**x²,y²,z²はこの像に入らない。**

## 3. Split SFとconstant-fixed-root C2の全空性

既constant-root必要条件を全split flagsへ適用した後に残るのは

* split000、h=xy、p=q=dlogz。
* split001、h=unused z、p=dlogy,q=dlogx。

だけ。残りの全split flagsは既active-component proofで空。
最初のprofileはX=Y=(-1,1,0)、κ=0。
q(C)=0はC_z=0、curl kernelはC∥∂zだからC=0、全空。
二番目はX=(1,0,−1)、Y=(0,−1,1)、κ=−1で、
`2κ+4Y−X=(-3,-6,3)`の全成分が非零。やはり全空。

C2のfixed rootがconstant uの場合、G3=u(Au²+2Buv+Cv²)。
Cはquadratic、R支持、q=dlog(f/C)。C=Li²なら三residues activeで空、従ってC=LiLj。
branch pairをx,yとしT=AC−B²のnonsquare partはxy。
degree四より全Tはconstant倍のxyz²、x³y、xy³。

C=xyなら全TがCで割れ、Bもxyで割れる。
constant normal shearでB=0とするとAはz²,x²,y²のどれか。
これらu³ coefficientsは§2のactualD3像に入らず空。
C=xzではmod zでT=x³y又はxy³はsquareになれないため除かれ、T=xyz²だけ。
Bもxzで割れるのでshear後A=yz。
C=yzはbase swapで同様。
唯一remaining profileは

`G3=z u(yu²+xv²)`, `p=(1/2)(dlogx+dlogy), q=dlogy`。

X=(1/2,-1/2,0)、Y=(1,0,-1)、κ=1/2。
curl diagonalは(9/2,3/2,−3)でkernel零だからconstant-root sourceは全空。
これでC2 constant-fixed-rootの全common-factor/zero complementsも処理した。

## 4. C2 moving fixed root degree一：F2とprimitive性

nontrivial C2でfixed root degreeは零又は一だけ（degree二ならresidual quadraticがconstantでsplit）。
degree一ではG3=ell F、Fはlinear binary quadratic、disc Fは二R linesの積。
normal GL2とbase component choicesによりF=xu²+yv²と取れる。
T=(F=0)⊂P2×P1はsmooth F2で、base projectionはquadric-cone double coverのminimal resolution。
H=C0+2f、normal P1のclassはf、div(ell|T)のclassはC0+3f。
vertex p=(x=y=0)上の全exceptional C0を保持する。

Rのprime preimagesはC0、二branch ruling fibers Fx,Fy、及びZ=(z=0)のirreducible section classH。
resultantがR支持なので、ell divisorはこれらの非負整数和。
classC0+3fとなる全場合は

* C0+b Fx+c Fy、b+c=3。
* Z+Fx又はZ+Fy。

前者でC0を含むことは、ellの二linear coefficientsがpで零、すなわちx,yだけに依存することを意味する。
parametrizationx=−v²r,y=u²rではell=r*g3(u,v)、g3は任意binary cubic。
support Fx,Fyだけならg3はsingle monomial。
`g3=βu³+δu²v−αuv²−γv³`という四coefficient表示から、その場合ellのrowは一成分だけで、common linear factorを持つ。
primitive degree一に反する。
後者の唯一sectionはz*u又はz*vで、これもnonprimitive。
従って**全C2 moving-fixed-root sourceは空**。
third line zがvertex pを避けるというtriangleの仮定をここで使っており、concurrentには適用しない。

## 5. C3 branch全三line：全conic pencilsとsquare image矛盾

Σ:w³=xyzはnormal cubic 3A2 surface。
三root mapsのdegree budget、Hodge/parity、RRによるconic-pencil→actual lineの全証明は
`PLANE_CONIC_LINE_M3_ROOT_COVERS_AND_EXCLUSION.md` §§6–7と同じであり、
ここでもH=−K、H²=3、χ=1、Du Valだけなので全仮定を満たす。
root mapsはdegree二conic pencilsで、Σ上の全linesからのhyperplane pencilに尽くされる。

lineをbase lineMへprojectionすると、w|Mはlinearでxyz|Mがcubeでなければならない。
Mがcomponent x/y/zならramification lineが各一つ。
他Mでは三coordinate-line intersectionsが一つに重なる必要があるが、nonconcurrentなので不可能。
従って全linesはこの三本だけ。
例えばw,zのpencilのnormをzで割ったleadingは

`G3=xy u³+z²v³`。

normal GL2で変えた全leadingも同じprofileを与える。
βはdiagonalなのでz²というpure cubic coefficientは§2のactualD3像に入らず、全sourceは空。

## 6. C3 branch二line：F3上の全root maps

branch pairをx,yとする。normal compactificationはP(1,1,3)、minimal resolutionはF3。
H=C0+3f、H²=3、HC0=0、Hf=1。
三root mapsはdeck orbitなので各H-degree≤2である。
F3上でfixed curvesを除いたmovable classD=aC0+b fはnef、b≥3a、HD=b≤2。
従ってa=0、D=f又は2f。
全root mapsはr=s/tを通るrational functions φ(r)、degreeは一又は二。
deckはr→ζr（ζ³=1、ζ≠1）。

degree一はnormal PGL2でφ=rへ直り、leadingは
`G3=h(xu³+yv³)`、hはRのline component。
degree二ではC(r)/C(φ)の唯一involutionτ∈PGL2を取る。
U上rootsは相異なるのでφ(r)=φ(ζr)はr∈C*で不可能。
従ってτ(r)=ζrの全solutionsは0,∞のみ。
tracefree representativeτ=[[a,b],[c,−a]]を使うとsolution polynomialは
`ζc r²−a(1+ζ)r−b`。
全rootsが0,∞ならmonomialであり、invertibilityからb=c=0、τ=-rだけ。
よってφはnormal PGL2を除きr²で、leadingはx²u³+y²v³。

このsquare coefficientsは§2のactualD3像に入らず空。
degree一でもh=x又はyならsquare coefficientが現れ空。
**h=unused third line zだけが残る**。
その全profileは

```
G3=z(A x u³+B y v³), A B≠0,
p=(dlogx+2dlogy)/3,
q=(2dlogx+dlogy)/3.
```

全leading amplitudes A,BをP1 sourceへ残す。
normal diagonal stabilizerでnormalizationすることは可能だが、固定P0上の全parallel spaceは二dimである。
このwholeprofileのhigher-source G上界は未証明。

## 7. Generic21の全flagsをcurlで排除する

02/10は既constant-root active-component proofで全空。
00ではh=Li²が空、h=xyだけ残るのでp=dlogz。
q=(q1,q2,q3)、sum q=1とすると

```
2κ+4Y−X=(3(1−2q3),−3(1−2q3),6(q1−q2)),
−2κ+4X−Y=(−3(q1+1),3(q2+1),−3(q1−q2)).
```

u-rootの非零bad Cを得るにはq_i=0かつ最初のdiagonal entryが零。
従ってqは(0,1/2,1/2)の三permutationsだけ。
しかしv-rootではp(C')=0がC'_3=0を要求し、第二diagonalの第一/二成分はこの三値で零にならない。
C'=0でDにv因子。よって全00は空。

01ではδ=L、h=H。
既residue比較によりH=Lならρ=−1又は2、H≠Lならρ=0又は2だけである。
actual character pairsとcurlを全て書く。

* L=z,H=z：p=(1,1,−1)、q=(0,3,−2)又は(3,0,−2)。
  最初のqではcurl=(30,6,0)、q-zero axis xのcurl非零。後者もx/y交換で同じ。
* L=z,H=x：p=(0,1,0)、q=(0,2,−1)又は(2,0,−1)。
  curlsは(9,−6,−9)、(−3,−18,3)。それぞれ唯一q-zero axisのcurlが非零。

従ってC=0となり、全01も空。
全marked component choicesはbase S3 permutationsにより上の全casesへ入り、他periodをgeneric値として捨てていない。
これで**全generic21 triangle sourceは空**。

## 8. Triple root：全periodから四diagonal profilesへ

G3=h u³、hは全quadratic monomials in x,y,z。
q=2p−dlog(f/h)とし、exponents h_iを用いr_i=1−h_i、sum r=1とする。
qのnonzero bad Cにはq_i=0が必要。
そのaxis上でcurl entryは

`(3/2)(3−h_i)(q_j−q_k)`。

3−h_iは1,2,3のどれかで非零なのでq_j=q_k=1/2。
従ってqは(0,1/2,1/2)の三permutationsへ全periodから落ちる。

off-diagonalr=α_Zのflatnessはχ(Z)=0、χ=p+dlogh。
χ_i=(q_i+1+h_i)/2は全て非零なのでZはdiagonal tracefree field。
rはclosed logarithmic、p−qもclosed horizontal非零。
triangleの二closed formsのwedgeは非零なのでflatnessはr=const*(p−q)を与え、
constant normal unipotentでrを消せる。
nonlog extensionを初めから零としたのではなく、全field equationとperiod条件から除いた。

q=(0,1/2,1/2)とし、base permutationsを除く全四pairsは

| h | p | q |
|---|---|---|
| x² | (−1/2,3/4,3/4) | (0,1/2,1/2) |
| y² | (1/2,−1/4,3/4) | (0,1/2,1/2) |
| yz | (1/2,1/4,1/4) | (0,1/2,1/2) |
| xy | (0,1/4,3/4) | (0,1/2,1/2) |

全leading targetは各h u³の全scalar lineである。
これらのhigher源を、fixed amplitude一点で代用しない。

## 9. 全五first symbolsと次の最弱finite target

全SFは§6のC3 unused-line一profileだけ、全generic21は空、全tripleは§8の四profiles。
全first symbolsはdiagonalで、rawX,Y及びnormalκuvは§2の全式で具体化出来る。
全bracketsはtoric logcanonical、従ってraw strict Poisson。
D2零のκを戻しているのでprojective rank≤2のactual first symbolsである。

固定plane上のtriangle choices六、anchor amplitude一、distinct normal charactersのGL2 conjugacy二から
firstsymbol source dimension upperは全五ともB=9。
全30 P1と9 P2、fullG3 target spaces（SFは二dim、tripleは各一dim）を保持し、
higher fiber dimension upper d≤10が得られればplane六とtensor scalar一を戻してB+d+5≤24となる。
このfinite全数学coverが確認されるまで、新しいoperator/minor計算は始めていない。
全concurrent三線、非被約R、m二/m四及び全E帰属は本稿の結論外である。
