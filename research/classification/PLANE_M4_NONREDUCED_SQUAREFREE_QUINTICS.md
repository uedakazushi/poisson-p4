# 非被約R・m四squarefree leadingの全実quintic有限化

2026-11-07、主義務G。producer全source幾何。新六full Casimir源の次元・rankは有限jobで判定する。
前のgeneric211空性とは別のsquarefree補集合を対象とする。

## 1. 非特別jはRのmultiplicityと無関係にimage-rank一

rank二plane S⊂SingD、actual Dは被約既約quintic、位数四とし、
G4∈H0(O(1))⊗Sym4 C2のgeneric binary quarticがsquarefreeとする。
parallel tensor G4/fからclassical invariantsの比J²/I³はconstantである。
I,Jが両方非零の場合、homogeneous degrees二/三とUFDから

`I=aL², J=bL³, Δ=cL⁶, a b c≠0`

となる。Lはnonzero base linear form。
parallel discrの零性からLはSupp Rの一成分だが、Rが3Lである必要は無い。
R=2L+Mでもこの証明を使う。

重要なのは、connectionのpole supportではなくG4自身のdiscriminantがLだけに支えられることである。
したがって**P2\L=A2全域**でroot incidenceがfinite étale degree四となる。
もとの接続がM上に高位極を持っていても、G4のroot coverはM上も正則étaleであり、
接続の延長を仮定せずこの大きなaffine openを使える。

A2のsimply connected性とPic零から、四rootはprimitive polynomial vectorsで表され、
各pair determinant及びoverall scalarはunit、すなわちnonzero constantとなる。
G4のbase polynomial degreeは高々一なので、四root-factor degreesの和も高々一。
少なくとも三rootはconstantであり、そのうち独立な二rootとのdeterminantがconstantなら残るrootもconstant。
従って全G4はaffine chartでconstantで、homogenizationからG4=L*g4。

よって全actual D=L*g4+g5はL,u,vの三変数だけに依存する。
既約性によりexact essential-threeであり、既受理の全八23成分に入る。
これはsimple-pole又はlog化を使わない、全Rに有効なweak recognitionである。

## 2. 特別jの全stable pencilと非被約Rの支持

I=0又はJ=0の全linear quartic spacesは、
`PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_PENCILS.md`§§2–3のpure algebraでimage-rank二以下となる。
この二節のWitt/Hankel/全symmetric-pencil証明はR被約性又はpole-orderを使わない。
rank一は前節と同じessential-three認識へ送る。
I=J=0はgeneric squarefree条件に反する。

rank二の全stable pencilsはI0の三相異J-root型、I0の二重+単純型、
J0のsecant型とcompression型で尽くされる。
I0三相異型のdiscriminantは三used linesを持つため、非被約cubicのsupportに入らない。
他三型は二相異used linesを必要とするのでR=3Lにも入らない。
従ってR=3Lの全squarefree m四源はrank一であり、既八23成分へ尽くされる。

R=2L+Mで残る全leading formsは、independent base forms A,Bとnormal u,vを用いて

`I0 mixed: A u⁴+Buv³`、

`J0 secant: A u⁴+Bv⁴`、

`J0 compression: A u³v+Buv³`

の三つである。全nonzero amplitudesはA,Bの選択に含まれる。
二lineのどちらが二重か、unused components及び全normal inverseを保持する。
以下ではR又は接続を固定せず、これらの全actual Dを直接扱う。

## 3. 全binary quintic項を戻した六actual正常形

いずれもD=G4+g5(u,v)である。
A,Bへnormal linear formsを加えるactual linear shearsは、g5の該当四monomialsを任意に吸収する。

**I0 mixed。** u5,u4v,u2v3,uv4を消すと残る項はαu3v2+κv5。
κ零ならDがuで割れるので非零。
normal diagonal、A/Bの補償scaleと全体scaleでκ一、α零/非零を0/1へ戻せる。
全二正常形は

`F_I,ε=A u⁴+Buv³+εu³v²+v⁵, ε=0/1`。

**J0 compression。** u4v,u3v2,u2v3,uv4を消した後はγu5+κv5。
γ零ならv、κ零ならuがDを割るので両方非零。
nonzero diagonal scalesとfinite fifth-root choicesでγ=κ=1へ戻し、全一正常形

`F_C=A u³v+Buv³+u⁵+v⁵`

を得る。隠れたcoefficient ratioは残らない。

**J0 secant。** u5,u4v,uv4,v5を吸収すると、残る項はαu3v2+βu2v3。
両方零、一方だけ非零、両方非零の全三場合を保持する。
一方だけの場合の二orderingはu/v及びA/B交換で同じ。
二つのnonzero diagonal ratiosと全体scaleで、全三正常形は

`F_S,0=A u⁴+Bv⁴`、

`F_S,1=A u⁴+Bv⁴+u³v²`、

`F_S,2=A u⁴+Bv⁴+u³v²+u²v³`

となる。g5の両残係数零をD非既約と誤って落とさない。
以上の六形はAに関してprimitive linearなのでC上絶対既約、四gradientはconstant independentなのでessential-four。
actual逆はminus shear、nonzero diagonal inverseと有限swapであり、全higher Poisson sourceも同時に可逆に動く。

## 4. 必要最弱の全canonical sourceと有限target

各Fに対し全150quadratic係数を保持した
`LF={divΠ=0, Π#dF=0}`
を使う。canonical strict lift/Pfaffian–Casimir identityにより、actual D class[F]を持つ全projective Poisson sourceはP(LF)に入る。
vertex w dependenceを残し、full nonlinear JacobiやactualPf条件を勝手に解いたと主張しない。

Y1=−4A∂A+u∂u、Y2=−4B∂B+v∂vとする。
F_S,0はY1,Y2を持ちDfamily≤17、F_S,1は2Y1−3Y2を持ち≤18。
F_I,0は−4A∂A−B∂B+u∂uを持ち≤18。
他三形はstabilizerを仮定せずDfamily≤19を使う。
ここではessential-four base選択Gr(4,5)の四次元を既に戻している。

したがってLF dimensionsが順に8,7,6,6,7,6以下なら全像≤24に十分である
（順序はS0,S1,S2,C,I0,I1）。
あるいは完全LF上のactual Euler-Pfaffianが恒等零なら、source dimensionを小さくしなくてもrank四空性で十分。

一回のfinite jobで六つのQQ全1075×150 source、完全kernel/nonzero minor、
全source上のquadratic Pfaffian mapを構成する。
失敗ならwhole necessary sourceとfalse criterionを保存し、実rank四/成分geometryへ移る。
接続の全profile、長いjet、GB又はpoint farmへは戻らない。

## 5. 全m四へ戻す際の分岐

nonredRのsquarefree leadingは、rank一essential-threeと本六形で尽くされる。
generic repeated leadingでは22/4及びimage-rank一は同じ八源、31は全-Rの別二形、
211高rankは`PLANE_M4_TWO_LINE_COLLISION_EXCLUSION.md`でrank四空。
本六shapeのfinite source判定が成功し、各別作者reviewを受理して初めてnonredR全m四G非優越が従う。
reducedRの既存m四全packetと合わせた全R結論もその後に限る。
m二/三及び全rank二plane/Gは引き続き別義務。
