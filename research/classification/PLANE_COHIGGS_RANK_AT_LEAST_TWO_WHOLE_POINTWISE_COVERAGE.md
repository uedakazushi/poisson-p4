# Rank零planeのO(1) co-Higgs・rank M二以上：全点の既知被覆

2026-11-08、PRIMARY E、log_manuscript。producer全source証明の集約。
新component数の増分零。独立受理・global採録はrootの別義務である。
1054主稿／frontier、旧code／report／job、global／Gitは不変。

## 1. 定理の正確な射程

P4上の一般rank四Poisson構造で、actual退化五次Dが被約既約、
rank零plane Pを持つとする。
全multiplicityでconormal可換となる既local modular補題を用い、
traceless co-Higgs共通imageのsaturationがO(1)、
coefficient map M:C3→sl2のrankが二又は三なら、
全sourceは既catalogue成分の有限閉和へ入るか、指定Dの仮定で空である。

より正確にはrank M二の全sourceは指定既約quinticで空。
rank M三はm≥3がrank二以下、semisimple全m二がD可約、
nil m二は既rank二plane/m三の全pointwise閉和又は空へ入る。
全actual frame／normalizerの逆、全higherとparameter零境界を含む。
rank M一以下、O(0)、scalar、zero-symbol、他のE源へは外挿しない。

## 2. 全sourceの入口はどこで尽くされるか

[全first-jet入口](PLANE_RANK_ZERO_MULTIPLICITY_TWO_FIRST_JET.md)のco-Higgs分類と、
[全m独立review](PLANE_HIGH_M_RANK_AT_LEAST_TWO_SEPARATE_REVIEW.md)の可換conormal復活を使う。
rank M三ではM可逆なのでM(U)非零、semisimple／nilpotentの二whole source、全tau。
rank M二ではimage-planeのKilling直交lineがnonisotropic／isotropic、従って
W=span(F,E)又はBorel W=span(H,E)の二型。
M(p)零／非零を分けた旧全actual normalizerと全逆が完全なsourceを与える。
全30のP1、全九P2、Euler六商、元350式を保つ。

|全source|全点の結論|証明|
|---|---|---|
|rank M三、m≥3|全rank二以下|全tau零八核→六shear＋二decomposable方向|
|rank M三、m二semisimple|全D可約|shear-exact全tau＋半の全九P2 cubic因子|
|rank M三、m二nil、generic／負半|全D可約|旧全shear-exactとP0のc²因子|
|同tau一|既known union|全二parameterに第二rank二plane/m三|
|同tau負三半|全D可約|全一parameterでc²因子|
|同tau三半|既known union|全二parameterに第二plane/m三、全R境界|
|同tau半|既known union又は空|元lambda mu一行で二face、第二plane／全cubic因子|
|rank M二、M(p)零、m二|全空|全finite shear＋残二whole源の元square row|
|同m≥3|全空|旧全共通U／二次因子／三row矛盾|
|rank M二、M(p)非零、m≥3|全空|zero-Lax四源のrank八束＋非零tau全Borel w因子|
|同m二、Borel image|全空|全旧164patchのw-row零、唯一旧wallのc因子|
|同m二、nondegenerate image|全空|下記の全ND四family|

rank M三の全tau coverは
[全semisimple／nil一](PLANE_COHIGGS_FULLRANK_M2_SEMISIMPLE_AND_NIL_ONE.md)と
[残三whole壁](PLANE_COHIGGS_NIL_THREE_WHOLE_WALLS_COVERAGE.md)、
[高m全rank二](PLANE_COHIGGS_FULLRANK_HIGH_M_WHOLE_RANK_TWO.md)を用いる。
陽性帰属はactual第二Poisson planeのrank二／m三を全parameterで確認した後、
既受理[全m三pointwise theorem](PLANE_RANKTWO_M3_WHOLE_POINTWISE_COVERAGE.md)へ送る。
Dの式だけ又は次元上界から親を推論していない。

rank M二M(p)零は
[m二全空](PLANE_COHIGGS_KERNEL_M2_WHOLE_EXCLUSION.md)と既high-m独立review、
M(p)非零のhigh-mは
[zero-Lax四源](PLANE_COHIGGS_NONKERNEL_ZERO_LAX_WHOLE_RANK_TWO.md)と
[Borel非零tau全因子](PLANE_COHIGGS_BOREL_HIGH_M_INVARIANT_HYPERPLANE.md)を使う。
Borel m二全体は
[全four-family空性](PLANE_COHIGGS_NONKERNEL_BOREL_M2_WHOLE_EXCLUSION.md)で、
全inhomogeneous/free P2と全二parameter壁を含む。

## 3. nondegenerate imageのm二：四whole familyのうち三つは全shear

全M(p)非零の旧入口ではW=span(F,E)ならR=0で、
NDSS(s≠0)、NDSSN、NDNSS、NDNの四familyに尽くされる。
全normalizerと逆は
[旧nonkernel独立review](PLANE_RANK_ZERO_NONKERNEL_SEPARATE_REVIEW.md)を用いる。

NDSSのgeneric kernelは六、actual shear六、D2単射。
旧selected-minorの全有理根s=±1,±1/3、四二次field

~~~
s²=1/2, 1/3, −1/18, −1/72
~~~

でも同じ完全核／shear／単射性が別作者により確認されている。
各fieldの全複素embeddingを含み、一方の数値平方根だけを扱わない。
s零はm二の外であり、旧全high-m NDSS0へ戻して全rank二以下とした。

NDSSNとNDNSSも全kernel六＝actual shear六、D2単射。
[全finite-shear補題](PLANE_COHIGGS_SHEAR_EXACT_SOURCE_REDUCTION.md)で
全higherを元P0へ戻すと、全actual五次は非零定数を除き

~~~
q_M(z,w)*det(E_base,T,U).
~~~

q_Mはbase一次×normal二次の非零cubic、
detはm二条件から非零のbase quadraticである。
従ってこの全三familyは被約既約quinticにならない。
全P2を先に零に置いた議論ではなく、全finite変換後の単射性で零とした。

## 4. 最後のNDN全七kernelを一方向へ戻す

残るNDNはU=∂b、M=bF+cE、T=c∂a、K=−z²/2。
元P0の全非零括弧は

~~~
Paz=cz, Paw=cw, Pbz=cw, Pbw=bz, Pzw=−z²/2.
~~~

既全K1は七、actual shear六、D2単射。
旧Sの同じ六minorから全30×7のK=S A+N bを読む。
保存順kernel座標alpha0,…,alpha6に対し
t=−10alpha5/3+5alpha6/3で、全actual残方向は

~~~
Nab=bw, Nac=cw/5, Naz=−zw/5, Nbz=3w²/5.
~~~

全finite shear／global inverseと一緒に全P2を動かし、P1=tNへ戻す。
旧D2の同じ九inverseでFULL九P2の必要唯一候補は零。
この九式だけからPoisson性を推論せず、元全350式を全tで戻す。
元row64、quad(a,b,c,z)、monomial cw³の係数は正確に

~~~
−12t²/25.
~~~

したがって複素被約支持ではt零が必要で、P0のproper因子から全空。
より弱い確認として、全候補P0+tNのactual五次そのものも

~~~
D=c(−3t²w⁴/25+2tcw³/5+c²w²−bc z²)
~~~

とc因子を持つ。
全c-rowがc倍であり、選択式の後に残る元Jacobiを解かなくても
この全模型にH=(c=0) invariantがある。
全zero境界／full higher／normalizer逆を戻して唯一NDN sourceも空となる。

## 5. ONE最後のwhole source実行

[inline job](jobs/plane_cohiggs_ndn_existing_whole_source.json)を
90秒／1GiBでgate PASSし、元の固定NDN whole sourceを一回だけ確認した。
実測0.068081秒、最大RSS22,176KiB。
[report](../../reports/CLASSIFICATION_PLANE_COHIGGS_NDN_EXISTING_WHOLE_SOURCE_20261108.json)
には全K=S A+N b、FULL九Q=0、全350／actual126、
全六finite shearとglobal inverse、全150 roundtrip、実行sourceと入力hashを保存した。
新M1、rank、minor、point、profile、jet又はGBは作っていない。
標準Fraction疎多項式の既独立primitivesだけを再利用し、producer／CASはimportしない。

以上で表の全sourceが閉じる。
依存する新producer各proofの独立受理が揃った後にのみaggregateを採録する。
全co-Higgs型の低rank、他のplane源、全E、成分上界及び分類主定理はなお別の未解決課題である。
