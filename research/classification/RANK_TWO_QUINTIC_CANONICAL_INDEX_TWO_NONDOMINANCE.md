# 次数五・rank二・canonical index二正規化の全G非優越

2026-11-07、主義務G、producer plane_collision_review。
**純幾何producer証明、独立review待ち。** 全canonical index二枝の射影Poisson像の次元を24以下とする。
全点のE帰属、新成分、全分類の完成は主張しない。新CAS実行無し。
既公開の初回部分還元の証明を§1–6に再掲し、§7以降で全残余を閉じる。
未pushから復元したtriple/markedノートは依存に使わない。

## 1. 既受理の全model入口

πはP4のrank四Poisson、実Dは既約被約quintic、非正規支持Sはintegral非退化degree五、π|Sは一般rank二。
Dのnormalization ν:X→D、L=ν*OD(1)のsectional genusを一とする。
既独立受理のFukuma全birational還元のcanonical枝では

`X normal canonical Gorenstein, -KX=2L, L ample and globally generated, L³=5`。

これはterminal又はQ-factorialという仮定ではない。
length公式からδS1で他の二次元非正規支持は無い。
既受理のclosed inverse-residue一形式αとdualizing volumeのdivisorから、
導手Weil divisor C=−KX∼2L。
nodeならCは一prime又は二prime、cuspならC=2Γ。
node primeごとのαはsimple-log pole、cuspでは一般にdouble poleを残す。
rank零専用のπX vanishingをrank二へ移していない。

## 2. 有限quasi-étale coverの不存在とH1(Xreg,Z)

これは上記polarized Xだけに依存する補題である。
有限connected quasi-étale cover f:X'→Xのdegreeをmとする。
KX'=f*KXはCartierで、discrepancyの有限cover公式からX'もcanonical。
L'=f*Lはample/GG、-KX'=2L'、(L')³=5m。
basepoint-free Bertini/adjunctionにより一般B'∈|L'|はnormal Gorenstein canonical del Pezzo surface、
`-KB'=L'|B', KB'²=5m`。
そのDu Val最小解消はsmooth weak del Pezzoで、有理性と最小有理曲面のK²からKB'²≤9。
従って5m≤9、m1であり、非自明な有限quasi-étale coverは無い。

特にClXの非零finite torsionはcyclic coverを作るので存在しない。
U=Xregとし、preimage(SingX)がSNC divisor Eとなるprojective resolution Z→Xを取る。
canonical Fanoのrational singularityとvanishingからH1(Z,C)=0。
全exceptional componentのclassesはH2(Z,C)で独立である。
実際、class零の実linear combinationはrelative numerically trivial、両符号のnegativity lemmaで零。
SNC complementのGysin列からH1(U,C)=0となる。
H1(U,Z)は有限生成かつrank零。非零なら有限abelian quotientからUの非自明有限étale coverが出来、
complex Riemann existenceとXへのnormalizationで上のquasi-étale coverに反する。
よって

`H1(U,Z)=0, ClX has no torsion`。                       (1)

terminal del Pezzoの分類定理をcanonical Xへ使わずに得た結論である。

## 3. nodeの全closed formはlogarithmic pencil

single prime Cでは一般L²曲線上のresidue theoremが10λ=0を与え、
rank二restrictionによるλ非零に反する。non-split conductorもこの場合に含む。
従ってnodeではC=Γ1+Γ2、各Γi L²=5、residuesはλ,−λ（λ非零）。

U上でα/λからdivisor Γ1−Γ2のlocal dlogを引くと、O_U(Γ1−Γ2)のflat holomorphic connectionになる。
intersections等のcodim二で残る見かけのpoleはΩ1のHartogsで消える。
(1)によりmonodromyはtrivial、global nowhere-zero horizontal sectionとそのinverseを得る。
このsectionはnormal X上のreflexive sheafへcodim二Hartogsで延長し、
proper XのGAGAによりalgebraic。非proper UでのGAGAを仮定していない。
従ってΓ1∼Γ2。和2LとCl torsion無しから

`Γ1∼Γ2∼L`。

defining sections s1,s2∈H0(X,L)を取ると、α−λ dlog(s1/s2)はX上reflexive holomorphic一形式。
既受理のcanonical-Fano/GKKP消滅H0(ΩX^[1])=0により

`α=λ dlog(s1/s2)`。                                  (2)

## 4. cuspではdouble poleを保ち、exact pencilになる

C=2Γから(1)を用いΓ∼L。
完成環K[[t²,t³]]⊂K[[t]]で、SのPoisson invarianceはmixed bracket{a,t}をtで割る。
base二変数のbracketはrank二なので一般に零でない。
volumeがt^-2を持つためαはnormal方向double pole、tangent方向simple poleを持ち得る。
閉性を展開すると

`α=d(a/t)+λ dt/t+regular`

で、λはconstant。唯一primeのpoleに一般L²曲線のresidue theoremを適用して5λ0、従ってλ0。
(1)とdivisor ΓのGysin列を使う。classΓ=LはU内のcompact一般L²曲線とdegree五で交わるので非零、
ゆえにH1(U\Γ,C)=0。αの全periodは零でglobal primitive fを得る。
local上の表示からfはΓにsimple pole以下で、O_X(Γ)のreflexive Hartogs延長とproper GAGAによりrational。
従ってΓ=(s0=0)として

`f=s1/s0, α=d(s1/s0)`、s0,s1∈H0(X,L)。                (3)

rank二からs1はs0に比例しない。
一primeだからdouble poleまで消えると誤って主張していない。

## 5. 全normalization moduleのweighted symmetric三行列

canonical vanishingとSerre dualityによりH1(X,mL)=H2(X,mL)=0、全整数m。
一般B∈|L|はcanonical Gorenstein del Pezzo degree五でh0(B,-KB)=6。
H1(OX)=0とrestriction列からh0(X,L)=7。
F=ν_*OXはintegral hypersurface D上のtorsion-free rank一ACM sheafである。
finite dualityは

`Hom_D(F,OD(2))≅F`

を与える。導手のgenerator qによるmultiplication β(a,b)=q abはsymmetric invertible formである。

新たに直接読解した一次資料：A. Beauville, *Determinantal hypersurfaces*, Michigan Math. J.48(2000),39–64,
<https://math.univ-cotedazur.fr/~beauvill/pubs/det.pdf>、Theorem B、印刷44–45頁（PDF index5–6）。
同定理はsmooth hypersurfaceに限定せず、integral hypersurface上のtorsion-free ACM sheafと
invertible symmetric formを仮定し、homogeneous entriesのsymmetric resolutionを与える。
本件はd5,t−3で正確にその範囲に入る。

最小resolutionを0→F1→F0→F→0とする。dualityからF1=F0*(-3)。
F0のdegree零generatorはunit一つ、degree一generatorは7−5=2個。
degree a≥2のgeneratorがあれば、dual F1にdegree3−a≤1のrelationが現れる。
最小性からそのcolumnはunit行へのlinear relation以外零であり、
ODのannihilatorはdegree五なのでそのlinear relationも零。injective resolutionに反する。
従って

```
F0=O⊕O(-1)^2, F1=O(-3)⊕O(-2)^2,
M=[[C3,A2,B2],[A2,l1,l2],[B2,l2,l3]], detM=nonzero_constant*D.
```

unitを第一generatorに取る。symmetric resolutionのadjugate pairingから

`q=l1*l3-l2²`

は上で使った導手generator（非零constant倍）であり、ν*qのdivisorはCである。
conductor idealはF(-2)で、そのquotientはCM pure surface。graded depth列も同じなのでSはarithmetically CMである。
δS1でgenerically reducedなので、そのconductor schemeはintegral Sそのもの。
従ってqはSを含むquadratic。非退化integral Sのためrank≤2は不可能、よってqはrank三。
特にl1,l2,l3は独立である。

なおこの同定はterminal-onlyのKuznetsov–Prokhorov分類に依存しない。
同論文Definition1.1/Theorem1.2(v)はterminal限定であることを直接確認し、本証明へ流用していない。

## 6. Vertex lineと全scheme構造

V=(l1=l2=l3=0)⊂P4はline。V上ではMのlower blockは零、detM0なのでV⊂D。
各p∈Vでfinite normalization algebraのfiberはnonzero unital algebra、unit1は零になれない。
Mの第二/三relationはA2(p)*1=B2(p)*1=0を与えるのでA2=B2=0、第一からC3=0。
従って全M entriesはI_Vに入り、M|V0。

T=ν^-1(V)をscheme-theoretic preimageとする。finite base changeから

`ν_*OT=OV⊕OV(-1)^2`。                                (4)

T→Vはfinite flat degree三、h0(OT)=1、h1(OT)=0、pa(T)=0。
resolutionの全global sectionsによりR_Xはdegree零unitとdegree一二generatorでR_D-moduleとして生成される。
ゆえにsection ringはdegree一生成、全|L|はXをP6へ閉埋込みする。
(4)からh0(T,L)=4でrestriction H0(X,L)→H0(T,L)はsurjective、kernelはspan(l1,l2,l3)。
従ってTはP6内P3をscheme-theoretically spanするarithmetically CM degree三genus零curveである。
reducible/nonreducedである可能性を残し、generic twisted cubicだけで尽くしていない。

## 7. 行列が保持する全extension dataと実Poisson像

以下ではlower blockを`[[u0,u1],[u1,u2]]`、上のentriesをC,A,Bと書く。
target P4の座標は固定しておく。従ってAは**全**I_V(2)、Bも全I_V(2)、Cは全I_V(3)から始まり、
vertex方向の係数やnormal方向のextensionを零にしていない。
Beauvilleのself-dual resolutionは、既に選んだF0→Fのsurjectionからそのdual comparisonを作る証明である。
したがってunitと選んだ二つのmissing degree一generatorを保持してsymmetric Mを取れる。
任意のmarkを含めたい場合、まずF0のdegree一basisをそのmarkに変え、その後relation側をdual basisにする。
unit以外のdegree一generatorはunitにtarget linear formを足す自由も含む。

β(a,b)=q abのscalarを固定し、adjugate pairingをこのβに一致させる。
qの非零scalarを変える場合はMの共通scalarとmarkの尺度を調整出来る。必要な平方根の選択は有限個である。
以下の恒等式はすべてnormalization algebra上の等式である。

```
adj(M)00=q=u0*u2-u1²,
adj(M)01=B*u1-A*u2=qY,
adj(M)02=A*u1-B*u0=qZ,
adj(M)11=C*u2-B²=qY²,
adj(M)12=A*B-C*u1=qYZ.
```

Y,ZはF0の二つのdegree一generatorを表す。
これらはS5=C[x0,...,x4]上の**module generators**であり、free basisではない。
1,Y,Zをfree basisとして係数比較する操作は使わない。
ただしdegree四以下のtarget polynomialがD上で零なら、Dのdegreeが五なのでliteral polynomialとして零である。
後のdegree四恒等式はこの理由で正当化される。

行列とmarksが定めるDとclosed formの組から、ambientのprojective Poisson tensorは高々一つである。
Dはanticanonical hypersurfaceなので、residueによるdualizing volume Ω_Dは非零scalarまで一意。
Dのsmooth dense openで`α=i_(π|D)Ω_D`はπ|Dを一意に決める。
restrictionのkernelは

`H0(P4, Λ²T_P4(-5))=H0(P4, Ω²_P4)=0`

であるから、[D,α]を固定すると[π]も一意である。
nodeの非零residue λ又はvolumeのscalarは[π]では追加parameterにならない。
Y/Z又はZ/Yは上のcofactorからtarget rational functionとして得られるので、
この対応は有限個のcoefficient chart上でalgebraicである。
具体的にはDとcofactorでdenominatorを払ったrestrictionの線形方程式をπの150係数に課し、
projective Euler同値を取り、唯一解が存在するchartで解けばよい。
したがって行列parameter族の次元上界を実Poisson像へ移せる。
以後、実際にはnormalizationやPoissonを与えない行列もsupersourceに含めてよいが、
像のfiberを引く際は**実際のadmissible source全体**で変換が保たれることを示す。

## 8. 低次数のambient Casimir補題と全従属marksの排除

**補題。** 既約quintic Dを実退化因子とするgeneric rank四πについて、
非定数`g=P2/ell²`（P2 quadratic、ell nonzero linear）がπ|DのCasimirになることは無い。

`η=ell*dP2-2P2*dell`はEuler contraction零のdegree二係数一形式、すなわちΩ¹_P4(3)のsection。
π#ηはT_P4(3)のsectionで、仮定からD\(ell=0)上零、従って全Dで零である。
よってπ#η∈H0(T_P4(-2))=0（Euler sequence）。
πがdense openでsymplecticなのでη=0、gは定数となり矛盾する。
これはD上のCasimirをambient rank四へ持ち上げる際の次数を正確に使った補題である。

node single-primeとnon-splitの場合は§3で排除済み。
split nodeではq=s1*s2、α=λdlog(s1/s2)。
各siがtarget linear system W=ν*H0(P4,O1)に入れば、Γiの像Sがhyperplaneに入って矛盾する。
mod Wで二markが従属なら

`s2=c*s1+ell`, c≠0, ell∈W nonzero。

ell=0なら二primeが同じになるので無い。
`Y'=s1+ell/(2c)`と置くと

`Y'²=q/c+ell²/(4c²)=P2`。

s1/s2がCasimirなのでY'/ellもCasimir、従ってP2/ell²がambient型のCasimirとなる。
q rank三よりこの関数は非定数、補題に反する。

cuspではq=s0²、α=d(s1/s0)で、s0∈Wは同じ非退化性で不可能。
s1=c*s0+ell、ell∈Wならell=0はα=0でrank二に反し、ell≠0なら
s0/ellがCasimirとなる。その平方q/ell²も非定数ambient型Casimirで矛盾する。
従ってnode/cuspとも、二markはmod Wで独立の場合だけが残る。
特にこれら二markをmissing generator Y,Zに選べる。
これはT reduced/nonreducedやXのterminal性によるcase切捨てではなく、全従属sourceの排除である。

## 9. Node：固定lower blockの全十一parameter source

Y=s1、Z=s2を取りq=YZと尺度を合わせる。
adjugateのdegree四恒等式から

`A*B-C*u1=q²`。                                      (5)

u1で剰余を取る。u0,u2は独立なので、polynomial UFD
C[u0,u2,v,w]において

`Abar*Bbar=u0²*u2²`。

両者はhomogeneous degree二であるから、全場合は有限なi=0,1,2に対して

```
A=alpha*u0^i*u2^(2-i)+u1*LA,
B=alpha^(-1)*u0^(2-i)*u2^i+u1*LB,
C=(A*B-q²)/u1,
alpha∈C*, LA,LB∈H0(P4,O1).
```

右辺の除算はpolynomialの割切れであり、u1≠0というtarget開だけに制限していない。
C∈I_V(3)も(I_V²:u1)=I_Vから従う。
従って全固定lower block sourceはdimension≤1+5+5=**11**。
LA/LBの全零係数、iの両極端、特殊因子化、nonreduced Tを含む。
normal性やπのJacobiが追加条件を課せばdimensionは減るだけである。

## 10. Cusp：全sourceと一つのcubic scroll

Y=s0、Z=s1を取りq=Y²とする。
degree四恒等式は

`C*u2-B²=q²`。                                      (6)

u2で剰余を取るとBbar²=−u1⁴。
従って二つのsign ε=±1をすべて保ち

```
B=ε*i*u1²+u2*L,                L∈H0(P4,O1),
C=(B²+q²)/u2,
Astar=A-ε*i*u0*u1-u1*L∈I_V(2),
Z0=Z+L.
```

u2の除算もpolynomialとして成立し、target上の除外chartを作らない。
全Cはこの式で決まり、Lは全五係数、AstarはI_V(2)の全十二係数を保持する。
determinantは直接展開で

`D=u0*q²-2ε*i*u1*q*Astar-u2*Astar²`                  (7)

となる。LはDから消えるが、primitive `(Z0-L)/Y`には残るので五parameterを落としてはいけない。

normalization X⊂P6上で

```
f1=Y²+u1²-u0*u2=0,
f2=u1*(Y+ε*i*u1)+u2*Z0=0,
f3=Astar+u0*(Y+ε*i*u1)+u1*Z0=0.
```

最初はq=Y²、残りはMのgenerator relationsである。
(Y−ε*i*u1)f2をf1で簡約すると
`u2*((Y−ε*i*u1)*Z0+u0*u1)=0`。
Xはintegralでu2非零なので
`f4=(Y−ε*i*u1)*Z0+u0*u1=0`。
f1,f2,f4はmatrix

```
[[u2,       Y-ε*i*u1, -u1],
 [Y+ε*i*u1, u0,        Z0]]
```

の三つの2×2 minorsである。
五つのactive座標(u0,u1,u2,Y,Z0)内では、これはSegre P1×P2のrank二hyperplane section、
すなわちsmooth cubic scroll S(1,2)である。
実際、六matrix座標の間の唯一linear relationは
`m21-m12+2ε*i*m13=0`で、対応する2×3 bilinear formはrank二。
残る二target vertex座標を加えると、P6内normal ACM degree三fourfold

`W4=S(0,0,1,2)`、vertex=P1

となる。全XはW4とquadric f3の交わりに入る。
f3はW4上非零である。三つのminorのconstant linear combinationにはmonomial u0*Yが無く、f3にはあるからである。
従ってW4∩(f3=0)はdegree六、pure dimension三で、degree五prime Xのほかに
degree一のeffective residual divisor R=P3を持つ。multiplicity X≥2はdegreeから不可能。

## 11. Residual P3の全補集合、全rのrank七

W4内の任意のP3はvertex P1を含む。
さもなければvertexからのlinear projectionの像はP2又はP3となり、nonplane scroll S(1,2)内に入れない。
従ってP3はscroll上のlineのconeである。
smooth S(1,2)=F1のlineは、ruling fiber f又は唯一negative directrix C0の二型で尽くされる。

directrix型を排除する。
二つの一般target-independent hyperplanesでP6をP4へ切り、vertexとSing Xを避ける。
W4の切断はF1、H=C0+2f。
Xの切断はsmooth integral degree五curve、adjunctionでK=(KX+2L)|=0、genus一である。
残余lineがdirectrix C0ならこのcurveはclass2H−C0=C0+4fで算術genus零となり矛盾する。
従ってRはruling型だけである。

rulingを(s:t)∈P1で書き、そのP3のfiber座標をU,V,v,wとすると

```
u2=s²*U,                 u1=-s*V,
u0=t²*U+2ε*i*t*V,        Y=s*t*U+ε*i*s*V,
Z0=t*V.
```

まず**例外ruling s=0**を処理する。
u1=u2=Y=0でu0,Z0,v,wは独立なので、f3|R=0はAstar∈(u1,u2)を強制する。
全導手idealは

```
I_S=(q,g1,g2),
g1=-ε*i*u1*q-u2*Astar,
g2=u1*Astar-q*L.
```

このときplane (u1=u2=0)上で三generatorがすべて零。
integral degree五surface Sが同次元planeを含むことは不可能である。
従ってこのrulingの全sourceは空であり、parameterを反転して黙って捨てていない。

残りはs=1、t=r∈A1の**全r**で

```
u2=U, u1=-V, u0=r²*U+2ε*i*r*V,
Y=r*U+ε*i*V, Z0=r*V,
f3|R=Astar(u0,u1,u2,v,w)+r³*U²+2ε*i*r²*U*V-r*V².
```

I_V(2)（dimension十二）のrestriction像は正確に

`span(U²,U*V,V²,U*v,U*w,V*v,V*w)`、dimension七。

u2=U,u1=−Vだから、このrank七はrの値に依存しない。
右辺のforcingはこの像に入るので、各rのaffine解fiberはdimension五である。
従って「あるruling Rでf3が零」というAstarの全必要sourceはdimension≤1+5=**6**。
r=0を含み、有限exceptional rootsやgeneric rank仮定は無い。
normal Xへのextensionが何であっても、全十二係数Astarから導いた必要条件なのでこの上界に含まれる。

## 12. 全target parameterと真のpresentation fiber

target P4内のprojective rank三quadricsのdimensionは

`dim Gr(2,5)+dim P(Sym² C³)=6+5=11`。

これは全target座標での値であり、後からPGL5の24を加えてはならない。
有限affine chartsでqのscalar representativeを一つ固定する。
固定qについてq=u0u2−u1²となるordered lower triplesは三次元orthogonal-isometry torsor。
従ってその全frame freedomはdimension**3**であり、qの共通scalarを重複して数えていない。
complementary vertex座標v,wはrestrictionを説明するためのlocal frameに過ぎず、
Astar,LA,LB,Lは最初から固定targetの全polynomialなので追加parameterではない。

nodeは§9よりpresentation source≤11+3+11=**25**。
各admissible presentationに対し

```
Y'=tY, Z'=t^(-1)Z, t∈C*,
M'=T^(-transpose) M T^(-1), T=diag(1,t,t^(-1)),
u0'=u0/t², u1'=u1, u2'=t²*u2,
A'=A/t, B'=tB, C'=C.
```

この変換はq、det M、normalizationのunit/algebra、node marksのdivisorsを保持する。
`dlog(Y'/Z')=dlog(Y/Z)`なので実際のprojective π像も同じ。
u0,u2非零よりstabilizerは有限、全fiberに一次元orbitがある。
これはsupersourceだけの対称性ではなく、全admissible sourceをpreserveするgenerator changeである。
従ってnodeのprojective Poisson像≤25−1=**24**。

cuspは§11よりpresentation source≤q11+lower frame3+Astar6+L5=**25**。
ここでは全admissible source上のadditive generator change

```
Y'=Y, Z'=Z+bY, b∈C,
u0'=u0-2b*u1+b²*u2, u1'=u1-b*u2, u2'=u2,
A'=A-bB, B'=B, C'=C
```

を用いる。これもdet-one congruenceでq、det M、normalization algebraを保持し、
`d(Z'/Y')=d(Z/Y)`を保つ。
二つのsign εの形へ戻してもB'=ε*i*u1'²+u2'*L'という同じ全sourceに入る。
u2非零なのでstabilizerは自明で、全fiberに一次元Ga orbitがある。
従ってcuspのprojective Poisson像≤25−1=**24**。

残余rulingを補助parameterとして数えたincidenceはpresentation sourceへ射影してdimension≤25を与える。
fiberを引くのはその射影像のうちのadmissible presentationsであり、特定ruling chartの非不変部分ではない。
したがってGaがrを動かす場合にも、不変な全sourceの最大dimensionから一を引く論理で有効である。

## 13. 全結論と残る分類義務

以上で、normal canonical Gorenstein index二normalizationを持つ全degree五rank二枝は
射影Poisson像dimension≤24となり、未知G成分の一般元（dimension25–27）から排除される。
single-prime/non-split node、dependent marks、両cusp signs、residual directrix、例外rulingを全て処理した。
nonterminal X、非normal Γ、nonreduced T、vanishing extension coefficientsは一切除外していない。
terminal-only del Pezzo分類、generic三次曲線の議論、formal jet、有限体実験は使っていない。

これは全点の既知成分帰属を証明する命題ではなく、E被覆と全分類の完成は引き続き別義務である。
新規性監査と独立reviewはproducer完成時点では未完了。
job gate及びhash manifestは保存整合性のみを確認し、数学的reviewの代用ではない。

