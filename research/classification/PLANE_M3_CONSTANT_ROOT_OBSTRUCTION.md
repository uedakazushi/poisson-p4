# Plane位数三・constant normal rootの全actualD障害

2026-11-07、PRIMARY G、producer plane_collision_review。
**新しい純証明、独立review前。** 新matrix、minor、GB、jet又はCASは無い。
全P1三十/P2九を保持した一つのactualD coefficientだけから全源の必要条件を得る。
reduced conic+line及びtriangleに対する以下の空性・有限character還元が新結論である。
concurrentとnontrivial SF coversの全完了は主張しない。

## 1. 仮定とconstant-rootの全gauge

actualDは被約既約quintic、S=P2はrestriction rank二、normal位数は三。
fはπSのcubicで、P_fのbase bracketを
`{a,b}=f_c,{b,c}=f_a,{c,a}=f_b`とする。projective chart c=1のanchorは3f∂a∧∂b。
G3のconstant normal rootをuとし、他normal coordinateをvとする。
このrootはsimple/double/tripleのどれでもよく、R被約性も本補題には不要。

global tangent Euler gaugeとmixed tracefree gaugeを取ると、全P0は

```
P0_base=P_f,
P0^{i,u}=u X_i,
P0^{i,v}=u Z_i+v Y_i,
P0^{u,v}=k20*u²+k11*u*v.
```

X,Y,Zは全tracefree linear base fields。
constant rootの保存がu-columnのv項を零にする。
G2=0とu-columnのこの形からnormal bracketのv²項は3f倍でD2へ現れるので零。
off-diagonal Zとk20,k11は捨てていない。

全higher source P1,P2を戻す。u-columnのうちuで割れない項は正確に

`P1^{i,u}|_(u=0)=v² C_i`, Cはconstant base vector。

他の全P1を固定せず、P2にはnormal columnが無いことだけを使う。
C=0ならfull u-rowの全成分がu倍となる。
hyperplane u=0でprojective tensorは三次元tangent space内に入りrank≤2だから、Dはuで割れる。
従って本対象では**C≠0が必要**。

## 2. 普通のnormal characterとの符号とfactor三

covector connectionのtriangular matrixをβ=[[p,0],[r,q]]とする。
constant normal uの保存からこの形で、qはcomplement diagonal characterである。
mixed field Yとの全関係は

`q=(1/3)dlogf−(1/3)α_Y`, `α_Y=i_E i_Y(da∧db∧dc)/f`。   (1)

例えばc=1で

```
{a,v}=v(Y_a−aY_c+f_b)+u(...),
{b,v}=v(Y_b−bY_c−f_a)+u(...),
q=(f_a−Y_b+bY_c)da/(3f)+(f_b+Y_a−aY_c)db/(3f),
α_Y=((Y_b−bY_c)da+(aY_c−Y_a)db)/f.
```

これでintrinsic affine quotientのEuler補正を含めて(1)が得られる。
raw mixed bracketsだけからqを推測していない。

## 3. 全actualDのv三乗係数

actualDのnormalizationを`E_total∧P∧P/2`の係数とする。
u=0でu-columnはv²C、v-columnの最低項はvY。
normal v³係数へ寄与するのは

* Eulerのv∂vとP_f∧v²C∧∂uからのC(f)。
* base Eulerとv²C∧∂u、vY∧∂vからの−det(E_base,C,Y)。

だけである。P1 base-baseはnormal degreeを更に一上げ、他のP1 normal columnとの積もv四乗以上。
P2はnormal columnを持たず、P0のZ/Kはuを含む。
従って全sourceのliteral coefficientは

`[D3]_(v³)=C(f)−det(E_base,C,Y)=3f q(C)`。             (2)

最後の等号は`det(E_base,C,Y)=(i_E i_Y vol)(C)`と(1)から従う。
normal root uによりG3|u0=0なので

**q(C)=0、C≠0**。                                     (3)

これが全constant-root obstructionである。
Jacobiのchosen sliceやfinite primesは使わず、full higher sourceのactualD coefficientを比較した。

## 4. R被約の場合の全active components

R=product f_iがreducedならqはclosed simple-pole character、
`q=sum ρ_i dlog f_i`, `sum deg(f_i)ρ_i=1`。
(3)をprime f_iで比較すると、ρ_i≠0ならf_i|C(f_i)。
右辺のdegreeはf_iより一小さいので

`ρ_i≠0 ⇒ C(f_i)=0`。                                 (4)

irreducible conic Qはnondegenerateなので、C(Q)=0ならC=0。
従ってconic+line R=LQでは、唯一許されるcharacterは

`q=dlogL`, `C(L)=0`。                                 (5)

triangle R=L1L2L3では三lineがindependentだから、三residuesが全て非零ならC=0。
従って少なくとも一ρ_iが零でなければならない。
concurrentでは全lineを殺す共通translation方向が存在するため、この論法を空性へ拡張しない。

## 5. Repeated21の全nonconcurrent結果

`PLANE_REDUCED_CUBIC_M3_ROOT_GEOMETRY.md`の全flagsとcharactersを用いる。
21のmoving degree1 characterηはR\L支持、Euler一。

### conic+lineの二型

00ではh=Q又はL²、p=dlog(f/h)、qは全Euler一character。
uとvの両方がconstant rootsなので、(5)をそれぞれcomplement qとpへ適用する。
h=L²はp=dlogQ−dlogLにQ-residueがあり空。
h=Qでは必要条件はp=q=dlogLである。

01ではh=L、δ=L、η=(1/2)dlogQ。
constant double-root uのcomplementは
`q=(3/2)dlogQ−2dlogL`、Q-residue非零なので全空。

02ではB=Qならq=3dlogL−dlogQ。
B=L²なら前noteの全normal-derivative還元からη=dlogQで、q=2dlogQ−3dlogL。
両方Q-residue非零、全空。

10ではconstant simple root vのcomplementは
`p=dlog(f/L²)=dlogQ−dlogL`、Q-residue非零で全空。

よってconic+lineの**全21 moving flagsは空**。
残00はh=Q、scalar ordinary connection β=dlogL Iだけに必要還元される。
この最後のfull higher sourceの上界又は帰属はまだ結論しない。

### nonconcurrent triangle

00でh=Li²ならp=−dlogLi+dlogLj+dlogLkの三residuesが非零なので空。
h=LiLjだけが残り、p=dlogLk、qの少なくとも一residueが零。
全qはこの三つのwhole affine period linesの和として保持する。

02でB=LiLjならq=3dlogLk−dlogLi−dlogLj、全residues非零。
B=Li²ならη=n dlogLj+(2−n)dlogLk、n=0,1,2なので
`q=−3dlogLi+(1+n)dlogLj+(3−n)dlogLk`、やはり全非零。
従って02は全空。
10ではconstant simple-rootのcomplementがdlog(f/L²)、三residues非零で全空。

01でδ=L、h=Hとする。残二lineをM,Nとしη=ρdlogM+(1−ρ)dlogN。
H=Lなら
`q=−2dlogL+(1+ρ)dlogM+(2−ρ)dlogN`。
よってρ=−1又は2だけ。
H=Mなら
`q=−dlogL+ρdlogM+(2−ρ)dlogN`。
よってρ=0又は2だけ。H=NはM/N交換で同様。
bは全linear formを保持し、零residueでCの方向が増える補集合も落とさない。
ここではperiodが有限値へ落ちたことだけを主張し、full P0/higherを一つの代表点で置換しない。

## 6. Fully split SFの全nonconcurrent結果

split SFの四degree flagsには少なくとも一constant normal rootがある。
000ではβ=dlog(f/h) I。
conicではh=L²が空でh=Qのみ、triangleではh=Li²が空でh=LiLjのみ残る。

001ではconstant u,vに対するp=dlog(f/(ha)),q=dlog(f/(hb))を両方検査する。
conicは相異なるline componentsを持たず存在しない。
triangleではh=a又はbなら一方のcharacterの三residuesが全非零で空。
したがってhはunused third lineでなければならず、p=dlogb,q=dlogaへ一意に落ちる。

002ではconicの{A,B}={Q,L²}の一方にQ-active characterがあり空。
triangleでquadraticが平方Li²なら対応characterの三residuesが全非零。
両方reduced quadraticsなら三components内でcoprimeに出来ない。
従って全002は空。

011でb1=b2=Lならcomplement q=dlog(f/L²)で、conicはQ-active、triangleは全三activeとなり空。
b1,b2 distinctの場合は既degree/resultant argumentでtriangleには存在せず、concurrentのみである。

従って全split SFに対し

* conic二型：000、h=Q、scalar β=dlogL Iだけ。
* triangle：000、h=two distinct linesのscalar型、及び001のh=unused thirdlineの一diagonal型だけ。

へ全源を必要還元した。SF quadratic/cyclic coversの残余をsplitと同一視しない。

## 7. Conic型SFの全C2 cover：fixed rootがconstantの場合

G3のrational fixed rootがconstant uで、残るquadratic root coverがnonsplitとする。
全形は`G3=u(Au²+2Buv+Cv²)`、A,B,Cはquadratics。
SFよりC及びT=AC−B²は非零、そのdivisor支持はRに含まれる。
ordinary complement characterは`q=dlog(f/C)`。
(5)からCはQのnonzero constant倍となり、尺度を合わせC=Qとする。
Tはdegree四でR=LQに支持されるから、全場合はconstant倍のL⁴,L²Q,Q²。
quadratic coverがnonsplitなのでTはC(S)で平方でなく、従ってT=τL²Q、τ≠0だけである。

`A=(B²+T)/Q`がpolynomialなのでQ|B²、prime QからB=bQ。
constant normal shearでBを零に戻すと、全leading sourceは

`G3=u(τ L²u²+Qv²)`, τ≠0。

全connectionは

`β=diag((1/2)dlogQ,dlogL)`

となる。normal diagonal scaleでτを調整出来るが、leading amplitudeをhigher sourceから消してはいけない。
これはnonsplit C2 coverの全constant-fixed-root branchであり、generic quadratic coverをsplitと扱っていない。
fixed rootがmoving degree一のC2 branch、又はconnected C3 branchは別残余として保持する。

## 8. Triple rootへの同じ入口

G3=h u³、β=[[p,0],[r,q]]、q=2p−dlog(f/h)。
constant root uに対して(3)を適用出来る。
conicではq=dlogLなので

* h=Q：p=q=dlogL。
* h=L²：p=(1/2)dlogQ、q=dlogL。

off-diagonalは全八dim sourceを保持してflatnessを課す。
r=α_Z（fixed nonzero scalar conventionはZへ吸収）、Z tracefree linearと書くと

`dr+(p−q) wedge r=0  ⇔  (p+dlogh)(Z)=0`。             (6)

実際dα_Z=Z(f) i_Evol/f²、Euler-horizontal ηについてη∧α_Z=−η(Z)i_Evol/f。
従って全flatness式から直接(6)が得られる。
これは全nonlog rを含むglobal linear-field条件である。

h=QではZ(f)=0となり、conic+lineのclosed horizontal form spaceは一dimなので
rはdlogQ−2dlogLのconstant倍。
h=L²では`4Q Z(L)+L Z(Q)=0`。
gcd(L,Q)=1からZ(L)=aL、Z(Q)=−4aQ。
Q nondegenerate、trace Z=0を取るとa=0、従ってZ(L)=Z(Q)=0。
やはりrは(1/2)dlogQ−dlogLのconstant倍で、constant normal unipotentにより消せる。
前者のnilpotent extensionはcoefficient零/非零の両方を保持する。

これでconic tripleの全first connectionsもscalar、scalar+nilpotent、diagonalの有限型へ必要還元される。
全higher源のG上界は依然別の義務である。

## 9. 残余と適用順序

新全空性はnonconcurrent reduced Rの上記moving/repeated・split SF branchesだけ。
concurrentではcommon translation Cを保持し、nontrivial SF covers、三線tripleのperiod/module、
最後のscalar/diagonal first symbolsの全higher liftは未解決として残す。
normal-root degrees、residue零の全補集合を復元した後に使うwhole coefficient obstructionであり、
P1のbad coefficient以外を零にするsliceは一切選んでいない。
このnoteを有限source probe又は全G完成の代用にしない。
