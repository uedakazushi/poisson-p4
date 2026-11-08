# Concurrent三線・位数三21/triple：全first-sourceと有限period壁

2026-11-08、PRIMARY G、producer plane_collision_review。純証明、独立review前。
本稿はfullhigherを計算せず、全first symbolを後続の一括判定に渡す。
全first-sourceの固定planeでのaffine次元は高々11。
全Poisson像のG≤24、E被覆、higher-source次元上界はまだ主張しない。

## 1. 入口、全marks、全gauge

actualDは被約既約quintic、rank二planeのnormal位数は三とする。
zero cubicは三相異なるconcurrent lines。
任意のmarked lineをL=a、残二本をM=b,N=a−bへ取る。
この表示は三本の全permutationsを許す。base第三座標cを保持する。
共通Pfはf=LMN、base bracketsは

`{a,b}=0, {a,c}=−a²+2ab, {b,c}=2ab−b²`。

generic leading cubicが21又はtripleならconstant normal rootが必ずある。
21の全primitive degree flagsは00,01,02,10、tripleのprimitive root degreeは零。
actual normal GL2でconstant rootをuとし、全raw first symbolを

```
{x,u}=uX, {x,v}=uZ+vY,
{u,v}=k u²+κuv,     x=(a,b,c)
```

と取る。X,Y,Zはtracefree linear base fields。
constant root obstructionから、fullP1のbad coefficient Cは非零でq(C)=0。
fullP1三十/fullP2九は未固定である。

普通のcovector connectionはβ=[[p,0],[r,q]]。
本稿の全p,qはc-independent logarithmic charactersなのでX,Yはc方向、κ=0。
以下のinverse formulaは全sourceで使う：

```
p=αL dlogL+αM dlogM+αN dlogN,   αL+αM+αN=1,
X=[(1−3αM)a+(1−3αL)b]∂c,
i_E i_Z(da∧db∧dc)=−3f r.
```

Yもpをqへ替えた同じ式。
E_total∧P0²=0のu²係数は正確に `3f(k−p(Z))` だから

`k=p(Z)`                                               (1)

がglobal quadratic symbolの必要十分なnormal K条件である。
K02=κ=0も全Schur係数から得る。
flat βと(1)のpolynomial性が揃えばU上でrank二horizontal Poissonとなり、
元の全polynomial projective Jacobiもdense openから恒等的に零。
raw liftのstrict Jacobi零は仮定していない。

## 2. Ore枝を先に完全に保存する

βがa,bだけを含むとZもc方向であり、k=0。
そのwhole P0は `∂c wedge W(a,b,u,v)` である。
以下を全てこのOre枝へ残す。

* 21の00：hは全六quadratic monomials in L,M,N、p=dlog(f/h)、qは全Euler一log character。
* 01/10でmoving rowの第三座標係数が零：全marked L,Hと全periodρ、全binary rowを保持。
* 02のsquare determinant L²：η=ρdlogM+(2−ρ)dlogN、全ρ。normal derivative条件から
  `A=M²+(ρ−2)LM mod L²`。normal shear/scaleで除いたL²とnonzero common amplitudeの逆もactual GL2である。
* 02のreduced determinant MNでAがbinary quadratic：§7中央族の部分。
* tripleでZがc方向：hは全六monomials、pは全二periods、q=2p−dlog(f/h)、
  `Z=(s a+t b)∂c` は全二係数。flatnessと(1)は自動。

01/10のbinary moving rowを一つの任意pointへ固定しない。
normal diagonal/shearを戻すか、元の全row係数を保持すればよい。
Oreという表示だけでfullP1/P2もc-independentであるとは主張しない。

## 3. Moving degree一の非Ore枝は共通Zρを持つ

moving rowを `b_row u+L v` と書く。
b_rowのc係数が零なら§2。非零ならactual c'=b_rowへ移り、anchor scalarを戻す。
逆はc=(c'−s a−t b)/d、d≠0で、全zero complementは§2に残る。
従って非Ore canonical rowを

`z=c u+L v`

と取れる。全nonlog offdiagonalは

`ηρ=ρdlogM+(1−ρ)dlogN, rρ=(cηρ−dc)/L`。

これは全ρ∈Cを保つ。inverse formulaから厳密に

```
Zρ=((1+ρ)a−3b)∂a+(ρ−2)b∂b+(1−2ρ)c∂c.             (2)
```

を得る。例えばi_E i_Zvolのda,db,dc係数は順に
`3(1−ρ)bc, 3(ρa−b)c, −3b(a−b)` で、−3frρと一致する。

### 3.1 flag01：H=Lは非Ore空、H=M/Nは一つの全ρ族

G3=H u²z、Hは全marked R line。
p=dlog(f/(HL))、q=p+ηρ−dlogL。

H=Lならp=dlog(MN/L)で

`p(Zρ)=ρ−2+3b/a`。

非constantであるため(1)を満たせず、この**全非Ore枝は空**。
H=Mなら

```
p=dlogN,
q=−dlogL+ρdlogM+(2−ρ)dlogN,
X=(a+b)∂c,
Y=((1−3ρ)a+4b)∂c,
k=ρ+1.                                               (F01)
```

H=NはM/N交換とρ→1−ρで含む。
leading targetは**全amplitudeの** `M u²z`。
元のb_row、anchor及びnormal frameを戻すactual変換を保持する。

### 3.2 flag10：normal swap後の一つの全ρ族

元のdouble rowはL old_u+c old_v、simple rootはold_v。
u=old_v,v=old_uへswapするとG3=u z²で

```
p=(1−ρ)dlogM+ρdlogN,
q=−dlogL+dlogM+dlogN,
X=((3ρ−2)a+b)∂c,
Y=(−2a+4b)∂c,
k=4ρ−2.                                              (F10)
```

Zは(2)。全ρでconstant kとなり、additional period cutは無い。
F01/F10のqはLにnonzero residueを持ち、他の少なくとも一本でもactiveなので
q(C)=0からbad Cは正確にC*∂cへ入る。C≠0も保持する。

## 4. F01/F10の全leading targetと有限period壁

u,zはordinary characters p,θz=q+dlogLを持つ。
full parallel polynomial G3をu,zで展開する。
flat direction ∂cにより係数はa,bだけに依存する。
z|L0=c uなので、異なるz-degreeのhighest c-powerは打ち消せず、
coefficientのL-poleを許すことは出来ない。
従ってregular termは `h_i u^(3−i) z^i`、i=0,1,2、deg h_i=2−iだけ。
各coefficientは

`dlog h_i=dlog(MN)−(2−i)p+(1−i)θz`                 (3)

を満たす。全nonnegative integral valuationsを調べる。

F01では

`h_i=M^[1+(1−i)ρ] N^[(1−i)(1−ρ)]`。

* 全ρ：M u²z。
* ρ=−1：さらにN²u³。
* ρ=0：さらにMN u³。
* ρ=1：さらにM²u³とu z²。

F10では

`h_i=M^[(2−i)ρ] N^[(2−i)(1−ρ)]`。

* 全ρ：u z²。
* ρ=0：さらにN²u³とN u²z。
* ρ=1：さらにM²u³とM u²z。
* ρ=1/2：さらにMN u³。

上記以外のperiodで追加regular parallel coefficientは無い。
F01ρ1は既SF square型の全三dim space、F10ρ1/2は既SF nonsplit型の全二dim space。
SFなのに一つのdouble targetだけへcutしない。各wallの全spaceを使う。
21/tripleに対象を限定する場合は、その実際のmultiplicity条件を全leading coefficientsへ戻す。

## 5. Triple全sourceのflatnessとglobal K

G3=h u³、h=L²又はLM（全marks）とする。
全βはp二periods、q=2p−dlog(f/h)、r全八dimから始まる。
Zをrのtracefree inverseとし

`χ=p+dlogh`, χ(E)=3。

flatnessは正確に

`χ(Z)=0`。                                             (4)

これはrをlogと仮定した条件ではない。
さらに(1)と(4)からp(Z)=−dlogh(Z)がconstantでなければならない。
prime residue比較により、Zはsuppχの全lineとsupphの全lineをそれぞれ保つ。
従ってunionの本数は一/二/三を尽くす。

三本ならZ(a)=s a,Z(b)=s bであり、(4)のEuler三から3s=0。
tracezeroからZcにc項も無く、§2の全Ore枝となる。
Zのa/b成分が零の場合も同じ枝へ入れる。

## 6. Triple・二line支持：同じZρと全resonance補集合

二本をM,N、omitted lineをLとする。
h=M^j N^(2−j)、j=0,1,2を全て保つ。
Z(M)=x M,Z(N)=y N、Zc=rM M+rN N−(x+y)c。
χ=λM dlogM+λN dlogN、λM+λN=3、(4)はλM x+λN y=0。

x=yならx=y=0でOre。非Oreではx≠y、normal diagonalでy−x=3へactualに直す。
ρ=x+2と置けばχ=(ρ+1)dlogM+(2−ρ)dlogN、全periodはCである。
canonical p,q,X,Y,kは

```
p=(ρ+1−j)dlogM+(j−ρ)dlogN,
q=−dlogL+(2ρ+1−j)dlogM+(j+1−2ρ)dlogN,
X=((3j−3ρ−2)a+b)∂c,
Y=((3j−6ρ−2)a+4b)∂c,
k=3j−2ρ−2.                                          (Tj)
```

ZはZρにrM M∂c+rN N∂cを加えたもの。
actual c-shift c→c+sM+tNはrM,rNを

`rM→rM+3(ρ−1)s, rN→rN+3ρt`

へ動かす。
ρ≠0,1なら両方を零へ戻せる。逆の全二shearsを失わない。
ρ0ではrNだけ、ρ1ではrMだけが残る。
nonzero残差はc-scaleとnormal diagonal、anchor scaleの組合せで1へ直せる。
従って全canonical補集合は

* 全ρ：Z=Zρ。
* ρ0：Z=Z0+N∂c。
* ρ1：Z=Z1+M∂c。

である。zero残差は第一族に含む。
one-line支持と重なるendpointも捨てず、§7との重複を許す。
これらのqもLと他少なくとも一本でactiveなのでbad Cは非零c方向。

### 全leading target

Z=Zρではu,zがdiagonal row frameなので、(3)から

`h_i=M^(j−iρ) N^[2−j−i(1−ρ)]`。

全ρでh u³、追加項は次だけ：

* i=1：ρ=j又はj−1で、対応するdegree一monomial h1 u²z。
* i=2：ρ=j/2でu z²。

全jのperiod wallsは−1,0,1/2,1,2だけ。
resonant nonzero残差では、u,z frameのoffdiagonalをdiagonal flat frameへ割り戻すと
nonzero constant倍のdlog(N/M)が残り、monodromyはnontrivial unipotent。
例えばρ0で追加R=N∂cなら

`L r_extra/N=(1/3)dlog(N/M)`。

ρ1でもMで割る同じ形（符号はframe規約）となる。
ρ0/1と整数jではdiagonal charactersのmonodromyはtrivialなので、
twistedSym³のunipotent invariantsはu³の一lineのみ。
従って両nonzero残差型のwhole leading spaceは正確にC*h u³。

## 7. Triple・一line支持と02 reduced determinantは同じ中央族

union一本ならh=L²、χ=3dlogL、p=dlogL、q=3dlogL−dlog(MN)。
(4)からZ(L)=0。
全tracefree linear field ZでZ(L)=0となるものは

`Z=3∇L cross ∇A`, A全quadratic modulo L²

と可逆に書ける。k=p(Z)=0。
全connectionは

`r=(2A dlogL−dA)/(MN)`。

これは02のreduced determinant B=MN、moving row m=A u+Bvの全式と同一。
primitive02にはgcd(A,B)=1を課すが、tripleには課さない。

qは三componentsでactive、bad C=λ∂c、λ≠0。
共通canonical transverse-offdiagonal補題
`PLANE_M3_CANONICAL_TRANSVERSE_OFFDIAGONAL_OBSTRUCTION.md`
SHA4fd66aea9bda0adc00c5b4fe500fbac21b142a109da1aed709c9dd0d9ce34347
を独立手確認済み。
strict uv²式のa/b projectionが(DZ C)_(a,b)=0を与え、Zb=−3A_cだから

`A_cc=0`。                                             (5)

全actual sourceは従って

`A=c ell(a,b)+Q2(a,b) mod L²`                         (J)

へ入る。全P1/P2を保持した必要条件であり、jet近似ではない。
ell0は§2のOre枝。
ell≠0ではactual c-shiftがQ2をell*linearで動かし、normal unipotentがAをA+tBへ動かす。
normal diagonalはAをscaleし、c-scaleはc ell部分だけをscaleする。

* ellがBのfactorでなければ、B mod ell≠0だから全Q2を除け、A=c ellへ戻す。
  ellの全P1 direction（Lも含む）を保持する。
* ellがBのfactorなら、Q2 mod ellはremaining lineの平方の倍数。
  zero/nonzeroを分けA=c ell+ε N_other²、ε0/1へ直す。
  primitive02にはε≠0、tripleではε0も保持。

全inverseは上記actual affine c-shifts、constant normal triangular変換とnonzero scales。
denominatorを持つ未定義PGLをzero branchへ適用しない。

### 中央族のwhole leading space

z=A u+Bv、B=MNとするとu,zのcharactersはdlogL,3dlogL。
diagonal flat frameu/L,z/L³でtwisted cubicを戻すと

`G3=α L²u³+βu²z+γu z²/L²+δz³/L⁴`。

LはBを割らず、z mod Lはvのnonzero coefficientを持つ。
Lのvaluationでδ,次にγが零となる。
よって**全A、全zero/nonprimitive境界で**regular parallel spaceは正確に

`span(L²u³, u²(Au+Bv))`、dimension二。

α,β両方を保持し、β0はtriple、β≠0ではprimitive性/共通factorによるdegree flagsを元に戻す。
係数AをL²で割ったquotientはconnectionだけの同定であり、失われたL² coefficientはこのαへ現れる。

## 8. 全first-source次元と、未了のwholehigher

concurrent reduced cubicのplane内familyはdimension五、anchor scalar一。
constant normal root uの選択はP1の一dimension。
固定uでのnormal complement changesは以下のr/A係数に既に含み、追加して数えない。

* 21/00：二normal root lines二＋simple character二、従って5+1+2+2=10。
* 21のbinary degree一又はsquare02：normalPGL2三＋period一、従って高々10。
* 非Ore F01/F10：root一＋period一＋全linear b_row三、従って5+1+1+1+3=11。
* triple Ore：root一＋p二＋Zc二、従って11。
* triple二line：root一＋χ一＋Zのeigen amplitude一/rc二、従って11。
  normalizing y−x=3で消費したamplitudeをこの数では明示的に戻した。
* 中央族：root一＋A mod L²かつA_cc0の四係数、従って11。

上界はzero/stabilizer jumpを許すsupersourceからの上界なので、
normalizerで固定したpointだけの計数ではない。
全marked選択は有限であり上界を変えない。

従って**全concurrent reduced-R m三repeated first-sourceのaffine dimensionは高々11**。
SFの先行finite-cover入口（まだ未列挙の候補を含む）もB≤9/10だから、
whole concurrent m三のfirst-sourceに同じ上界11が使える。
これはfullhigher sourceの次元ではない。

全P1三十/P2九を戻し、actualD3を各whole parallel targetへ合わせたfiberの上界d≤8が
一括に得られれば、plane六とprojective scalar一を戻して11+8+6−1=24が十分。
SFの小さいfirst-sourceでは既記録の緩いfiber上界でよい。
全350Jacobi/fullactualDへの実lift、period壁上の全kernel、wholeG上界は未実行・未証明。
本稿だけから新component又は既成分帰属を数えない。
