# Whole T37の単純component低次数leaf条件とGa9排除：別作者純査読

2026-11-08、`/root/generic_exhaustion`。**独立純読解・別手導出 PASS。**
新CAS、point、M1/operator、minor、profile又は局所変形jetは実行していない。
producerのimmutable V1とsimple-component scope補足を分離して受理する。
新component数／全分類完成は本receiptから加算しない。

## 1. 受理する正確なwhole条件

π0∈T37がrank四とし、actualD0のmultiplicity一のprime component B上で、
Bのsingular locus及び他componentを避けたp0を取る。π0(p0)rank二なら、
p0はdegree≤5のintegral Poisson invariant surfaceに属する。
**重複componentのreduced supportがsmoothなだけの点には適用しない。**
元 `T37_WHOLE_CLOSED_LOW_DEGREE_LEAF_SURFACE_THEOREM.md` 第4節の一般component表現は、
補足 `T37_LOW_DEGREE_LEAF_SURFACE_SIMPLE_COMPONENT_SCOPE.md` に従ってこの意味だけで採用する。
元V1／旧manifest bytesは変えない。

## 2. generic T37入力を逆推論せず固定する

canonical index二のadmissible source→T37という片方向定理から、
一般T37のnormalizationをcanonical Fanoと逆断定してはいない。
実際の固定canonical tensorは真completed T37 all-oneと全150/Eulerで同定済み。
同じall-oneをK4aのe一tensorへ移す別作者全150等式と、K4a非零全eが一actual24軌道で
その閉包がT37である独立定理を合成する。
したがってこの固定canonical tensorのPGL軌道自体がT37にdenseである。
必要なのはこのactual dense-orbit橋であり、normalizationのconverseは不要。

固定F,g1,g2に対するinverse-residue dlog(g1/g2)も、既root別Fractionが
全六二形式の元polynomial restriction identityとして確認した入力である。
そのためgeneric pointのπ_D annihilatorをrational pencilへ接続出来る。

## 3. HB baseとdegree boundの独立監査

q=ac−b²,A=c²+bd,B=a²+be、g1=bB−cA,g2=bA−aB、
F=Cq−cA²+2bAB−aB²、bC=AB−q²とする。
直接展開でbF=g1g2−q³。
I=(q,g1,g2)は[[A,B],[a,b],[b,c]]のmaximal minors。
qは既約、q0でも(a,b,c)=(0,0,1)でg1=-1なのでgrade≥2、determinantal boundでgrade2。
Hilbert–Burch resolutionは

```
0 -> O(-4)^2 -> O(-2) plus O(-3)^2 -> I ->0.
```

Hilbert numerator1-t²-2t³+2t⁴の二階微分は1で10、degreeS=5。
b0ではsupporta=c0のlineのみ。bunitではq0からa,cunit、
a g1+b g2=-AqによりI=(q,g1)、qでc、g1のe係数b²unitでeを消せる。
b1 chartはa∈C*,d自由だからsmooth integral。
HBのCM/unmixednessとgenericreducedからS全体はintegral、cycle multiplicity一である。

generic S上のcodimension二regular local ringでIはmaximal ideal、F∈I²。
κnonzeroの一般値についてGκ=g1−κg2は一次normal partを持つ。
従ってCI(F,Gκ)のS multiplicityは≥2。
より直接にはgeneric A,b,a unitでg1,g2がnormal parametersとなり、Fの最低項はunit*g1g2、
Gκで割るとunit*κg2²となる。弱い≥2だけで十分である。
Fは既約degree5、Gκは非零degree3なので全CI cycleはpuredegree15。
全Sを少なくとも2S差し引いた残effective cycleのdegreeは≤5。
S外の一般点pのκ=g1(p)/g2(p)を取ると、その点を通るintegral residual componentがある。
residueの全式によりこのcomponentはπ-invariantである。
よって一枚の例示surfaceではなく、genericDの一般点をdegree≤5 surfaceで覆う。

## 4. DVR/Hilbert specializationの独立監査

π0∈closureO37からgenericが上記denseorbitにあるDVR arcを取る。
D0非零なのでactual quintic familyはprimitiveなrelative degree5 Cartier hypersurfaceでflat。
単純Bの指定p0ではD0自身がsmooth、従ってD_Rはそこでrelative smooth。
完備局所環R[[u1,u2,u3]]で一般sectionを選び、generic p_Kをpencilのbase／bad loci外へ取れる。
必要なら有限field extensionを行う。これはmarked pointを別のspecial pointへ動かす操作ではない。

p_Kを通るintegral invariant surfaceZ_K一枚のHilbert polynomialを固定すれば、
projective Hilbert schemeのpropernessがflat Z_Rを与える。
一様な全degree≤5 Hilbert polynomial分類やChow tangency閉定理を仮定する必要はない。
p_R∈Z_RとZ_R⊂D_Rはgeneric inclusionとclosedness／flat torsionfreenessで延長する。

局所Hamiltonian derivationδ_Rとg∈I_ZRについてδgはgenericに零modI、
R-flat quotientのtorsionfreenessから全体で零modI。よってspecial idealもPoisson stable。
char零でstable idealのradicalもstable（g^nへn回微分しradical内で還元する）。
reduced ringのminimal primeへ局所化するとfieldなので、各minimal primeもderivation stable。

p0を含むreduced componentはrank二Hamiltonian distributionを含みdimension≥2。
fixed Hilbert polynomialのdegreeは二だから全special support dimension≤2。
従ってp0はtopdimensional invariant surface componentに乗る。
flat leading degree≤5がeffective cycleの各component degreeを抑える。
低次元embedded supportだけへp0が逃げる可能性もこのrank二論法で排除される。

## 5. Ga9 targetの全leafとdegree九を別導出

h1、F4=0でw=2y²−y⁴/4、三変数(x,y,z)のπは
(4−y²,−x+2yz+4w²,y)。標準volumeのkernelは

```
y dx+(x−2yz−4w²)dy+(4−y²)dz
 = d[xy+(4−y²)z−16y⁵/5+4y⁷/7−y⁹/36].
```

yと4−y²は共通零無し、従って全このopenでrank二。
Cを上のCasimir、T=(-x+yz)/4とすれば(y,T,C)は全polynomial coordinates。
全finiteC=κはintegral A²であり、一つのconnected leafである。

scroll resolutionP_P1(O4⊕O⊕O)、L tautological,F fibre、
(h,y,w)=λ(s⁴,s³t,2s²t²−t⁴/4),x=X,z=Zで、leaf closureは

```
s⁴tX+s³(4s²−t²)Z
 -λ[(16/5)s⁴t⁵-(4/7)s²t⁷+t⁹/36+κs⁹]=0.
```

classL+5F。三fiber-linear coefficientsにcommonbase zero無し：sunit,t0ではZ係数がunit、
tunit,sunitではX係数がunit、s0ではλ係数t9/36unit。
よって全finiteκでintegralな相対hyperplane、verticalsurfacecomponent無し。
L³4,L²F1、h1でbirationalなのでimageのdegreeはL²(L+5F)=9。
このmapはresolutionであり、λ0のvertexline contractionを持つので有限normalizationとは呼ばない。

もしGa9unit orbitがwholeT37内なら、F4のh1 smooth rank二点で§1のdegree≤5 surfaceが必要。
そのsurfaceは同点のA²leafを含み、同じdimension二でintegralだからleaf closureそのもの、degree9。
矛盾。したがって**Ga9unit全23軌道はwholeT37の外**。
その磁気零boundaryの全非所属やC_G∩T37=emptyを主張しない。

## 6. 同じwhole条件のDOUBLEへの適用も別導出一致

`CONSTANT_MODULAR_DOUBLE_CUBIC_NO_ALGEBRAIC_LEAF_SEPARATE_PROOF.md`を独立に読解し、
印刷されたcanonical全括弧から同じ三変数計算を得た。
h1の単純cubicでxunit,y=-(x+1)/x²、u=t-((x+1)(x+2)/x³)rとすれば、
全annihilatorはunit倍のdu+dx/x。
characteristic distributionは∂rとx∂x−∂u。
Laurent UFDでinvariant irreducible surface equationを取ると、∂r-invarianceでr非依存。
Darboux cofactorのx最大最小指数／u次数比較からcofactorconstant、
Σx^j p_j(u)の比較からp'_j=(j-c)p_j、従って唯一eigenpolynomialsはunits x^j。
非unit equationは存在しない。
よってこの単純cubicの一般rank二点にはinvariant surface自体が存在せず、
§1の必要条件と矛盾し、**DOUBLE forced orbitもwholeT37の外**。
DOUBLEのdoubleH componentへ§1を適用したのではない。

## 7. 保存と未主張

producer V1、simple scope addendum、全source/target旧bytesは不変。
本receiptは純別作者数学査読であり、CASや旧全matrixの再実行と偽らない。
T37以外の親排除、全prior一覧、新成分性／総成分数はrootが別依存と統合して判断する。
