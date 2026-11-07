# 次数五・rank二・canonical index二正規化：全pencilと三次曲線への還元

2026-11-07、主義務G、producer plane_collision_review。
**部分結果。** 以下はcanonical index二枝の全構造還元であり、全G≤24又は空性の完成ではない。
終点に残るreducible/nonreduced三次曲線を明示する。新CAS計算は行っていない。

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

## 6. vertex lineの全triple curveとintegral補集合の排除

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

cuspならq=s0²上constant倍。T reducedの場合、s0|Tの平方零からs0|T0。
よってs0∈span(li)⊂ν*H0(P4,O1)であり、Γの像Sはambient hyperplaneに入る。非退化性に反する。
従ってcuspはT nonreducedの場合だけに残る。

node splitならq=s1s2。T integralの場合どちらかsi|T0となり、同じ理由でSがhyperplaneに入り矛盾。
従ってnodeもT reducible又はnonreducedの場合だけに残る。

## 7. 厳密に残る数学

canonical index二全枝は、nodeではreducible/nonreduced、cuspではnonreducedの全minimal-degree cubic curve Tへ還元された。
その全normalization/projection/Poisson支持のdimensionをまだ押さえていない。
特異X、非normalなΓi、三次曲線のline/conic splitting、nilpotent thickeningは除外していない。
全sourceのG≤24又は新componentの不存在は未証明。

次の構造的入口は(4)の全cubic algebraと、rank三qのvertex lineに沿うPoisson moduleである。
新しいparameter点列、terminal分類の無断拡張、whole-source GB又はjet反復は行っていない。
本部分還元の独立reviewと新規性監査も未実施である。
