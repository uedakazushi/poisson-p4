# R二重直線＋直線・m三squarefree・connected三root coverの全first symbol

2026-11-07、主義務G。producer純幾何。current runner/jobsのfreeze中に保存する。
新operator計算又はG上界はまだ主張しない。

## 1. 対象と全root cover

actual Dは被約既約quintic、S⊂SingDはrank二plane、S沿い位数三、
πSの零cubicをR=2L+M、LとMは相異直線とする。
generic leading binary cubic G3∈H0(O(2))⊗Sym3 C2はsquarefreeと仮定する。
本noteはU=S\(L∪M)上のdegree三root incidenceが**connected**である全枝を対象とする。
root coverがsplit三成分又は一＋二成分、m三重根、m二/四は別補集合である。

Hamiltonian associated-gradedによりG3/fはfull conormal connectionのparallel tensor。
generic stable binary cubicのtwisted GL2 Lie stabilizerは零なので、このtensorが接続全体を一意に決める。
高位極及びnonlog性は最初から許す。

base coordinates(a,b,c)でL=(a=0)、M=(b=0)とし、chart b=1でt=a/b、r=c/bを取る。
U=Gm×A1、π1(U)=Zなので、connected degree三finite étale coverのmonodromyはcyclic三。
標準Kummer cover t=s³で全三rootが分かれる。
そのcoordinate ring C[s,s^-1,r]はUFDでPic零、unitはconstant*sの整数powerだけ。
各rootをprimitive polynomial vectorで表すと、pairwise determinantsとoverall coefficientはunitとなる。

各root vectorのr-degreeはunit倍に依らずwell-definedで、cyclic三actionにより全て同じd。
G3のr-degreeは高々二で、leading productの消滅はdomain性から起こらないので3d≤2。
従ってd=0、全G3はrに依存しない。これは一つのchosen cyclic exampleではなくconnected全coverの帰結である。

## 2. 全first connectionは一変数のdouble/simple-pole形

G3/fからの接続一意性によりordinary conormal connectionもtだけに依存し、dr成分は零。
global normal bundleはO(-1)⊗C²、Rのpole multiplicitiesはaで二、bで一である。
global H0Ω1(3)↔tracefree linear base fieldsの同型を使うと、
その全raw weight零representativeはEuler gaugeを除いて

`P0=∂c∧[a²∂a−2ab∂b+((aM1+bM2)n)·∂n]`, n=(z,w)

と書ける。M1,M2は初めは全2×2 constant matrices。
D2=0がnormal-normal quadraticの残三係数を一意に零へ戻す。
これは接続がc-independentだから得られる完全表示であり、一般の非connected枝へは拡張しない。

chart b=1では{t,r}=−3t²、{r,n}=[M2+t(M1+2I)]n。
従ってordinary conormal connectionは

`β=(C/t²+D/t)dt, C=M2/3, D=(M1+2I)/3`。

constant dt項はb=0で二位極を作るため、Rのb側一位極という正当なglobal degree条件に反する。
twisted binary-cubic discriminantのGL2 characterはdet²なので、parallel性から
`dlog Δ(G3/f)+2 trβ=0`。
左のrational unitのlog derivativeはsimple poleだけを持ち、従ってtrC=0。

## 3. 三algebraic Riccati rootsはC非零を許さない

全rootはalgebraic functionsであり、connected三monodromyのためnormal coordinateのinfinityが恒等的な一rootになることはない。
finite affine rootをq_i(t)、i=1,2,3とする。
接続のRiccati方程式を

`q'=A(t)+B(t)q+E(t)q²`

と書くと、任意pairについて

`(q_i−q_j)'/(q_i−q_j)=B+E(q_i+q_j)`。                  (1)

nonzero algebraic functionのlog derivativeはPuiseux座標で高々simple pole。
以下ではこのelementary事実だけを使い、irregular ODEの外部分類を仮定しない。

**Cがnonzero semisimpleの場合。** trC零なのでconstant normal GL2でdiagonalに出来る。
するとB=κ/t²+b1/t、κ非零、A=a1/t、E=e1/t。
e1零なら(1)のdouble poleが消えず矛盾。
e1非零なら全pair sumsのleading termは−κ/(e1t)、従って各rootのleading termはd/t、d=−κ/(2e1)。
Riccati右辺のt^-3係数はκd+e1d²=−κ²/(4e1)非零で、左辺のpole order二と矛盾する。

**Cがnonzero nilpotentの場合。** constant normal changeでRiccatiのdouble-pole項をA=κ/t²+a1/t、κ非零へ置く。
B=b1/t、E=e1/tとなる。
e1非零なら(1)から全pair sums、従って全q_iがPuiseux-holomorphicとなり、
Riccati右辺のκ/t²を左辺又は他項が消せず矛盾。
従ってe1零であり、方程式はlinearである。

b1零ならalgebraic性がa1零を強制し、全rootは−κ/t+constantsでmonodromyがtrivialとなる。
b1=−1ならnonzero κのlog項を生みalgebraic rootは無い。
他の場合はconstant normal translationでa1項を消し、全三rootは

`q_i=c0/t+d0 ζ^i t^b1`, c0 d0非零、ζ³=1

となる。connected三coverはb1のdenominatorが三、3b1整数かつb1非整数を意味する。
normal cubicはunit factor h(t)を用いて

`h(t)[(w−c0 z/t)³−d0³ t^(3b1) z³]`

である。全rootはC*上にpoleを持たないので、w³係数hはC*上nonzeroなpolynomial、h=constant*t^p、0≤p≤2。
z²w係数のregularityからp≥2、従ってp=2。
z³係数にはnonzero t^-1項が残り、消すには2+3b1=−1、b1=−1が必要。
これは上で除いたresonanceでもあり、connected三coverにも反する。
従ってC非零の全場合は不可能、M2=0となる。

## 4. 全実normal formsは三つ

C零ならβ=D dt/t。
connected三rootのfinite projective monodromyはnontrivial order三なので、Dのprojective nilpotent partは零。
constant normal GL2でM1=diag(λ,μ)へ直せる。
三rootはEuler solutionsのcyclic orbitであり、G3の全形は二extreme monomialsの和

`A0 a^p b^(2−p) z³+B0 a^q b^(2−q) w³, A0 B0≠0`。

p,qは0,1,2の相異なる整数。normal z/w交換でp<qとしてよい。
全parallel equationのcoefficientsは

`3p−4+2λ−μ=0`, `3q−4−λ+2μ=0`。

従って全first symbolは次の三actual正常形で尽くされる。

| (p,q) | (λ,μ) | 全parallel regular leading space |
|---|---|---|
| (0,1) | (3,2) | <b²z³, abw³> |
| (0,2) | (2,0) | <b²z³, a²w³> |
| (1,2) | (0,−1) | <abz³, a²w³> |

全sourceでは両leading amplitudesを自由に保ち、その積非零がconnected/squarefree条件である。
normal diagonalization、base ordered-line changes及びtensor scaleはactual invertibleで、
元の任意first symbolと全higher sourceへ逆に戻れる。
algebraic gaugeだけでPoisson全源をPGL同一視したのではなく、それを不可能性証明にだけ使い、
最後に実在するconstant normal changesで全representativesを得た。

## 5. finite jobの正確なtargetと残余

固定planeでordered lines四、πS nonzero amplitude一、normal eigenbasis二から、全first-symbol源dimensionは七以下。
全higher sourceはP1三十、P2九。
full weight一JacobiとD3が表の全二次元spaceに入る条件の完全核を求め、actual D3 imageを直接調べる。
imageがnonzero squarefree cubicを持たなければ全枝空。
持つ場合、fullweight二operator kernelを加えたhigher上界dが十二以下なら、
plane六とprojectivization一を戻して全像7+d+5≤24で十分。

一回のbounded三model jobで全350/30/9、全24D3係数、complete kernels/minorsと
leading-discriminant polynomialを保存する。
足りない場合は全sourceを保持し、元Jacobi又はactualD/rank四geometryへ移る。
新しいminor/jet/点列へ自動反復しない。
split cover、一＋二cover、m三repeated leading、m二及び他Rは未完了として保持する。
