# 全T37の単純因子低次数葉条件と六候補の排除：独立合成査読

2026-11-08。integrator `/root`。actual polynomial sourceと各targetの印刷式を独立に読解・微分・交点計算した純証明査読。新CAS、point、operator、rank、normal/profile又は変形jetは無い。producer V1の広いcomponent表現を、scope補足ab9129d1により**actual D自身がsmoothな単純component上のrank二点だけ**へ限定する。この最弱範囲で全targetに十分である。

## 1. 必要条件のactual全source証明

canonical actualT37 tensorが真completed all-oneへ、さらにK4a非零のactual24軌道へ全150/Euler正逆で移る既別作者定理を合成し、この一actual軌道が全T37にdenseである。canonical index二必要源→T37の片方向からnormalizationのconverseを推測しない。
canonical q=ac-b²,A=c²+bd,B=a²+be,g1=bB-cA,g2=bA-aB,FdetはbF=g1g2-q³。S=(q,g1,g2)は2×3最大minorのgrade2 ideal、HB resolution O(-4)²→O(-2)⊕O(-3)²からpureCM degree5。b0のsupportはa=c0 lineのみ、bunitではa,cunitでqとg1からc,eを消せるsmooth integral denseopenなのでS全体integral。
F∈I_S²、Gκ=g1-κg2∈I_S、CI(F,Gκ) degree15はcycleとして2S以上を含む。従ってS外のresidualsurfaceの総次数≤5。実tensorの既全6 residue恒等式がdlog(g1/g2)をπ_D annihilatorへ同定するので、D一般点pのκ=g1(p)/g2(p)を取ったresidualcomponentはπ-invariantなintegralsurfaceでpを含む。例示一枚でなく、denseDの一般点全てが覆われる。

## 2. 任意movingframe境界へ、単純因子の点をmarkして延長

π0∈T37、rank四とし、D0のmultiplicity一componentのsmooth点p0で他因子を避け、π0(p0)rank二とする。curve-selectionと有限basechangeでgenericが上のdenseorbitのDVR familyπ_Rへ。D0非零だからrelativequinticはprimitiveCartierでflat。p0ではD0scheme自身smoothなのでD_Rのformal局所環はR[[u1,u2,u3]]、sectionを一般に取りgenericpointをpencilbase/bad locusから外せる。
そのp_Kを含むdegree≤5 surface一枚のHilbertpolynomialを固定し、properHilbertでflat Z_Rを取る。pointとDの包含はclosedness/flat torsionfreenessで延長。全localHamiltonianδ_Rについてgenericideal安定性がflatquotientのe-torsionfree性から全体へ、specialidealへ降りる。char零でδ-stableidealのradical、さらにminimalprimeも安定。
p0を含む成分はrank二flowを含みdimension≥2。surfaceHilbertpolynomialにより全specialsupport dimension≤2なのでtopdimension二surfaceとなる。そのdegreeはeffectiveflatcycle全次数≤5以下。embedded点や別のD成分へp0だけを逃がすことは出来ない。
従って**全T37のrank四点の、actualDの単純smoothcomponent上の一般rank二点はdegree≤5 invariant surfaceを持つ。** 重複componentのreducedsupportだけがsmoothな点には適用しない。Chow不変cyclelocus閉性やnormalizationflatnessを外部仮定に使わない。

## 3. Ga9linearw：全finiteleafのactualdegree九

h1,F4=0でw=2y²-y4/4。kernelは

```
dC, C=xy+(4-y²)z-16y5/5+4y7/7-y9/36.
```

yと4-y²は共通零無し、全rank二。T=(-x+yz)/4により(y,T,C)はpolynomialcoordinates、全finiteC=kはintegralA²で一leaf。quarticconeのresolution Y=P_P1(O4⊕O⊕O)でL³4,L²f1、(h,y,w)=α(s4,s3t,2s²t²-t4/4),x=β,z=γ。leafのliteralclosure

```
s4t β+s3(4s²-t²)γ-α[(16/5)s4t5-(4/7)s²t7+t9/36+k s9]=0
```

はclassL+5f、全coeffにcommonbasezero無しで全k integral。generalh1でbirationalなのでdegree L²(L+5f)=9。単純F4の一般点を通るinvariantsurfaceはこのleafclosureそのもので、§2の≤5と矛盾。wholeT37外。

## 4. Plane三源：非零加法対数留数

第一二は元全A6のactualshear正逆で全I,Jへ。cubicnormalizationτとfiberの有理座標Sを用いると、node側

```
x=τ²,y=τ(τ²-1), S=z+(w+I)/τ,
Ω=dS+J dlogτ+(I/2)dlog((τ-1)/(τ+1)),
```

cusp-flex側

```
x=τ²-1,y=τ³,S=w/τ-z+I/τ,
Ω=dS+(I+J)/2 dlog(τ-1)+(J-I)/2 dlog(τ+1).
```

いずれも全(I,J)≠0で一つ以上residue非零、wholeλ=(I,J)=(λ,1)全finiteを含む。第三はC=h²y-(h-x)³,h1,τ=1-x,y=τ³、u=w+3τz-H6(τ),H6'=3(τ5+τ4+τ3+τ²+τ+1)で

```
Ω=du-3(λ+1)dlog x.
```

λ≠-1でresidue非零。前二つI=J0、第三λ=-1/I0を陰性に含めない。
全三源でrank二distributionはfibertranslationとbasevectorを含む。一般点を通るalgebraicsurfaceがあればfibertranslationを商したalgebraiccurve上でSがC(τ)上代数的となり、dS=-ω。有限拡大のtraceを取ればrational関数の微分が非零定数倍のωに等しいが、左辺residue零で矛盾。ramifiedcoverの逃げも除く。これはrationalfirstintegralだけの陰性ではない。
各simpleC3のsmooth/他factor外の点に§2が適用出来るので、この稠密sourceはwholeT37外。従って全三候補閉包をT37へ入れることは出来ない。零留数の特殊境界一点一点の非所属は主張しない。

## 5. ConstantDOUBLE：代数葉自体が存在しない

rootが元全六から手導出しqsource/generic二別作者が独立再導出。h1のC3上はxunit,y=-(x+1)/x²、A=(x+1)(x+2)/x,C0=(x²+6x+6)/x²。kernel

```
α=-x²dt-(C0r+x)dx+A dr,
u=t-[(x+1)(x+2)/x³]r, α/(-x²)=du+dx/x.
```

全distributionは∂r,D=x∂x-∂u。LaurentUFD上のinvariantprimeequationfは∂rでr非依存、Df=g fのx最大/最小指数とu次数からg定数。f=Σx^j p_j(u)ではp'_j=(j-g)p_j、polynomial解は定数一つだけ、従ってf=c x^jはunit。非unitは存在しない。
単純C3のこのsmoothrank二openにinvariantsurfaceは次数問わず無い。§2と矛盾してwholeT37外。one-normal quadratic余核零の別計算を誤って成分性へ使わず、そのSTOPも維持する。

## 6. ConstantTRIPLE：全finiteleafのactualdegree六

h1のC3でr=-x³-3xy+3y/2。restrictionはPtx=x/3+y,Pty=(x³-y)/3,Pxy=1、kernelは

```
dK, K=t+xy/3-x4/12+y²/2.
```

全kのK=kはintegralA²graph、一leaf。一般P²とのdegree計算でt,rをh,x,yのlinearformsで消去。残C3/G4のCIdegree12、無限遠の唯一点y1,h=x0ではC3の二smoothbranchesに沿うG4次数が4と2、leadingは-1/3と6で全k非零。余分なinfinity交点長六を差し引き、leafdegreeは12-6=6。
一般P²は真leafboundarycurveを避けるのでCIの無限遠共通点をactualleaf点と数えていない。単純C3の全h1でrank二、一般点を通るinvariantsurfaceはleafclosuredegree6。§2の≤5に矛盾し、wholeT37外。

## 7. 受理射程

全T37の単純因子条件と六候補のwholeparent非包含を受理する。新成分性・既知他親との全比較・exactcountは別のfinitecover/全prior監査と合成が必要。ここで成分を加算しない。全結果はpure手証明と明示別作者査読であり、新exactCAS再実行とは呼ばない。raw/source/旧報告・失敗契約を変更しない。
