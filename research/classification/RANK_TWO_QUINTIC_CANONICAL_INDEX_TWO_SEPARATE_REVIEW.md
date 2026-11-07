# canonical index二の全G≤24：別作者による全証明の監査

2026-11-07、PRIMARY G、独立review **PASS**。
対象は凍結producer `RANK_TWO_QUINTIC_CANONICAL_INDEX_TWO_NONDOMINANCE.md`（SHA0b3eeb0f…）と
その初回還元である。新CAS、producer import、数値point/prime検査は行っていない。
以下は全sourceの証明確認であり、全Poisson分類又はE membershipを完成したとの意味ではない。

## 入口とcanonicalの全補集合

仮定はintegral非退化degree五S、一般rank二restriction、既約被約quintic D、
canonical Gorenstein normalization Xと`-KX=2L,L³=5`である。
terminal、Q-factorial、smooth S、normal conductor component、reduced triple curveは仮定されていない。
Fukuma還元後のこの枝を入力とする範囲を保ち、rank零専用のvanishingを使っていない。

有限quasi-étale coverのpullbackではcanonical/Gorensteinとample GG index二が保たれる。
一般|L'|曲面はDu Val del Pezzoで、degree5m≤9からm=1。
例外divisorsのclassesの独立性はrelative numerical trivialityと両符号のnegativity lemmaから従う。
Gysinとrational-Fano vanishingによりH1(Xreg,C)=0、さらに有限abelian coverの不存在からH1(Xreg,Z)=0。
Cl torsionのcyclic coverも排除される。

single-prime/non-split nodeのresidue10λ、split nodeのresidue±λを一般L²曲線上で確認した。
flat line bundleのmonodromyをH1=0で消し、codim二Hartogs後の**proper X上**GAGAで
Γ1∼Γ2∼Lとα=λdlog(s1/s2)を得る。非proper smooth locus上のGAGAを用いていない。
cuspではnormal方向double poleを保持し、唯一primeのresidueを消した後、Gysinでperiodを消す。
primitiveのsimple poleをreflexiveに延長してα=d(s1/s0)を得る。double poleを誤って捨てていない。

## ACM、選択generatorを保つweighted symmetric matrix

全整数mのH1/H2 vanishingはKodaira-type vanishingと`KX=-2L`によるSerre dualityで両側を覆う。
h0(X,L)=7、F=ν*OXのACM性とrank一、finite dualityによるinvertible symmetric formを確認した。
Beauvilleの一次資料を独立に直接読解した：*Determinantal hypersurfaces*,
Michigan Math. J.48(2000),39–64、Theorem B、printed44–45（PDF index5–6）、
<https://math.univ-cotedazur.fr/~beauvill/pubs/det.pdf>。
同定理はintegral hypersurface上のtorsion-free ACM sheafとinvertible symmetric formを扱い、smooth限定でない。
証明終端の対称化は選択済みsurjection p:L0→Fを保持するため、unit及び実marksを固定したF0を使える。

degree≥2 generatorのdual relationはunit行の低次数annihilator以外に入り得ず、Dのdegree五に矛盾する。
従ってF0=O⊕O(-1)²、F1=O(-3)⊕O(-2)²、全entriesはC3,A2,B2と三linearformsで尽くされる。
adjugate pairingのq=lower determinantは導手generatorであり、CM purityとδ=1から導手schemeはintegral S。
q rank≤2ならSがhyperplane内に入るためrank三、三lower formsは独立となる。
vertex上のunitはzeroになれないので全M entriesがvertex idealへ入り、全extension coefficientsを保持する。
これらはterminal-only分類を使わない。

[D,α]はD smooth dense open上π|Dを決め、ambient restriction kernel
H0(Λ²TP4(-5))=H0(Ω²P4)=0によりprojective πも一意となる。
cofactorで分母を払った有限個のcoefficient chartからalgebraic partial mapとして扱える。
normalizationやJacobiを満たさないsupersourceを含めても上界は増すだけである。

## 全従属marksのCasimir排除

producerはg=P2/ell²を用いる。η=ell dP2−2P2 dellはΩ1(3)であり、
π#ηがD上消えるとH0(TP4(-2))=0からambient上消える。
一般rank四と非定数性に反する。degree五Hamiltonianに直接D因子を取る誤った次数論ではない。

nodeでmod W従属ならs2=c s1+ell、両marks非ambientからc≠0、ell≠0。
平方完成でY'²=q/c+ell²/(4c²)となり、Y'/ellの平方がこの禁止Casimirを与える。
cuspでs1=c s0+ellなら、ell=0はα=0、ell≠0はq/ell²が禁止Casimirとなる。
q rank三から非定数性が確かであり、全zero/特殊係数を含めて独立marksだけが残る。

## Nodeの全十一parameterと真C* fiber

独立marksをY,Zとしq=YZに合わせる。degree四identity AB−u1C=q²をu1で剰余すると
Abar Bbar=u0²u2²。polynomial UFDとdegree二から全i=0,1,2、α≠0、全LA,LB各五係数で尽くされる。
α=0は未処理chartでなくproduct identityにより不可能である。
u1による割算はpolynomialの割切れでtarget coordinateを反転していない。
従って固定lower sourceのdimension≤11、全vertex extensionも含む。

projective rank三qはkernel plane6＋quotient symmetric form5でdimension11。
qのaffine代表を固定するとlower ordered framesはO3 torsor、dimension3。
vertex補完座標はpolynomialsを記述するlocal frameに過ぎず、新たなsource parameterでない。
全presentation source≤25である。

Y→tY,Z→t^-1Zのdet-one congruenceはq,D,unit/algebra,αを保持し、
u0,u2独立からstabilizer有限なので**全admissible fiberに**一次元C* orbitがある。
従って実projective Poisson像≤24。supersourceだけの対称性を引いていない。

q scalar代表の正規化も確認した。q'=c qに対してa b=c、λ²=c³とし、
marksをaY,bZ、行列をλ diag(1,a,b)^(-transpose) M diag(1,a,b)^(-1)へ移す。
q'とmark productが同時に合い、Dと[α]はscalarしか変わらない。平方根は有限選択である。
したがってq scalarを追加で数える必要はない。

## Cuspの旧scroll上界も全範囲で成立する

q=Y²からCu2−B²=q²。u2剰余UFDは両signを保つB=εi u1²+u2Lを与え、L全五係数を残す。
Astar全I_V(2)十二係数を保ったdeterminant、三generator relations、u2をintegral X上で消して得る
第四relationを別に手展開して確認した。三minorの唯一linear relationはrank二なので
active surfaceはsmooth S(1,2)、全coneはnormal ACM S(0,0,1,2)である。
f3は非零quadricでdegree六cutを作り、degree五X以外のeffective残余にP3が必ず存在する。

W4内のP3はvertexを含む。含まなければprojection imageはP2又はP3でnonplane scrollに入り得ない。
scrollのlineは唯一directrixと全rulingだけである。
general codim二cutはvertex及びcanonical Xのdimension≤1のsingular locusを避け、
Xcutはadjunctionでsmooth genus一curveになる。
directrix残余ならclassC0+4fの算術genus零となり矛盾する。

exceptional ruling s=0は導手idealからSが同次元planeを含むこととなりintegralityに反する。
s=1,t=rの**全r**ではu2=U,u1=-Vで、全I_V(2) restrictionの像は
U²,UV,V²,Uv,Uw,Vv,Vwの正確な七次元空間。
forcingは常にその像に入り、affine fiber五＋r一からAstar source≤6。
r=0その他の特殊値にrank低下はなく、target openへの不当な制限もない。

q11＋frame3＋Astar6＋L5=25。
Z→Z+bYのdet-one congruenceはq,D,normalizationとd(Z/Y)を保持する。
u2非零からGa stabilizerは自明で、全admissible presentationの実π fiberに一次元orbitがある。
補助rulingが動いても、不変な全admissible presentation像でfiberを引くため上界24は正しい。
cuspのq scalar合わせもa²=c,b=1,λ=cで可能であり、[α]はscalarしか変わらない。

## 独立に得た短い強化と最終受理範囲

さらに全cusp sourceは空まで強められる。
別ノート `RANK_TWO_QUINTIC_CANONICAL_CUSP_EMPTY_SEPARATE_REVIEW.md` に全sourceの短証明を分離した。
導手idealのg1=−εiu1q−u2Astarから、integral非退化S上u2 generic非零なのでAstar|S=0。
qと独立quadricならdegree四CIにS degree五は入らず、Astar=cqならD=q²*linearで矛盾する。
従属marksの全排除と独立marksの全表示を合わせているので、normality/terminal/vertex curve等の追加開はない。
この強化を前提にしなくても、上の旧scroll/Gaによる全≤24は独立に確認済みである。

従ってproducerの全canonical index二G像≤24を**受理する**。
採録時はnode≤24、cusp空という簡約を使える。
全G一般元の排除であり、他degree・他normalization・node全点のE帰属や全分類完成は主張しない。
後続generic triple表、未push復元資料、旧kernel/jet反復は証明の前提にしていない。
