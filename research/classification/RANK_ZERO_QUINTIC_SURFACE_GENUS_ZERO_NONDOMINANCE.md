# degree五rank零支持・sectional genus零枝の全G非優越性

2026-11-07、主義務G。新producer幾何証明。独立review前はglobal完了主張へ混ぜない。

**命題。** rank四Poisson tensor πの実Dが既約被約quinticであり、
degree五integral surface S上でπが零とする。
Dの正規化(X,L)のsectional genusが零である分岐の全PGL5射影像はdimension≤24である。
Sのnormality・smoothness、plane残余上のπのrank、又はvertical divisorを仮定しない。
sectional genus一の分岐はこの命題の明示補集合である。

## 1. 正規化と唯一のplane残余

entranceの `Σδ_i degS_i=6−g` にdegS=5、g=0を入れる。
δ_S=1であり、他に一つのplane P、δ_P=1がある。二次元の非正規支持はSとPだけ。
plane PはPoisson invariantであり、そのrestriction rankは零又は二の双方を保持する。

`RANK_ZERO_SEXTIC_CONDUCTOR_EXCLUSION.md` §1–2と同じBertini/Höring Theorem1.2の
全normal-scopeにより、Xの完全L埋込みはdegree五minimal-degree threefold in P7となる。
Xがconeなら、そのlinear-projection image Dもconeで、degree五成分Sもconeとなる。
既受理の全高次数rank零cone排除とhyperplane内degree≥3補題で矛盾する。
従ってXはsmooth scroll

`X_113=P_P1(O(1)⊕O(1)⊕O(3))` 又は
`X_122=P_P1(O(1)⊕O(2)⊕O(2))`

だけである。quotient conventionでL=O_X(1)、fiber divisorをFとし、
`K_X=−3L+3F`、`C_X=C_S+C_P∼3L−3F` とする。
δ一の横断一般型はnode又はcuspで、C_S,C_Pのcoefficientはそれぞれ一又は二。
`L² C_S=10`、`L² C_P=2` である。

## 2. plane上の全degree二導手divisor

X上のeffective divisor A=aL+bFについて、一般P2 fiberへrestrictするとa≥0。
degree二なら `5a+b=2`。
a=0ならA∈|2F|。a≥1のglobal sectionは

`H⁰(P1,Sym^a(E)⊗O(2−5a))`

であり、Eのmaximal degreeは三以下なので、a≥2ならmaximal summand degreeは
`3a+2−5a=2−2a<0`、不可能。
a=1の場合、X_122ではmaximal degree2−3<0、やはり不可能。
X_113では唯一のsectionがあり、divisorは

`Q=P(O(1)⊕O(1))∼L−3F`

というsmooth quadric surfaceである。
従ってC_Pは全場合を含め `2F` class、又はX_113のQそのものだけ。
nodeなら前者は二つのfiber、cuspなら二倍の一つのfiberである。

## 3. Sのgeneric cuspはrational plane quinticのGauss次数に反する

Sがgeneric cuspならC_Sは偶数divisorである。
C_P class2Fでは `C_S∼3L−5F` は偶数でないので不可能。
残る可能性はX_113、`C_P=Q∼L−3F`、`C_S=2Γ_S∼2L` だけとなる。

しかし一般P2切断のplane quintic Cはrationalであり、Sとの五つの一般交点はordinary cuspとなる。
これは次のelementary Gauss boundに反する。

**Gauss補題。** rational integral plane quinticは五つ以上のordinary cuspを持たない。
normalization f:P1→C⊂P2をdegree五の三つのbinary formsで書く。
tangent-line mapの三Wronskianはdegree八である。
各ordinary cuspではdifferentialが一次で零となるため三Wronskianに一つのcommon linear factorを与える。
五cuspがあればdual curve C*のdegreeは高々三。
標数零でtangent-line mapはbirationalであり、bidualはCである。
これは座標で `g=f×f'` に対し `g×g'` が非零rational factorを除いてfとなることからも分かる。
C*のnormalizationもP1で、そのdegreeをe≤3とすると、同じWronskian構成で
bidualのdegreeは高々 `2e−2≤4`。
original Cのdegree五に矛盾する。補題が従う。

従ってこの全分岐でSのgeneric transverse typeはnodeであり、C_Sはreducedである。
Pがnode/cuspであることの両方は引き続き保持する。

## 4. 正規化bivectorのcodimension一でのregular性

πのD上のrational restrictionをXへpullbackしたπ_Xを使う。
S上でπが零であり、Sがgeneric nodeなので、前掲sextic noteのbranch quotient議論で
π_Xはreduced C_Sに沿って零となる。

Pのgeneric nodeでもbranch quotientによりπ_Xはregularである。
Pのgeneric cuspでは、normalization parameter tに対してlocal ringは
`A=K[[t²,t³]]⊂B=K[[t]]` である。
PがPoisson invariantなので、baseのregular rational function aに対し
`{a,t²}∈(t²,t³)`。
従って `{a,t}={a,t²}/(2t)` はBにregularで、むしろtで割れる。
base同士のbracketもAにregularである。
よってπ_XはP上でregularであり、π|Pのrank零/二を追加制限していない。
その他のcodimension一ではνは同型である。
smooth XでHartogsを使えば

`π_X∈H⁰(X,Λ²T_X(−C_S))≅H⁰(X,Ω¹_X(C_P))`。       (1)

global normalization-Poisson-extension定理を未確認で呼び出したわけではなく、
δ一の全node/cuspで必要なregular性を直接確認した。

## 5. 全bivector源の小ささ

relative cotangent完全列

`0→O_X(−2F)→Ω¹_X→Ω¹_X/P1→0`

を使う。
C_P=Q=L−3Fなら、relative部分をP2 fiberへrestrictして
`H⁰(Ω¹_P2(1))=0`。
base部分は `H⁰(O_X(L−5F))=H⁰(E(−5))=0`。
従って `H⁰(Ω¹_X(Q))=0` であり、この全caseは空。

C_P class2Fなら、relative部分のfiber H⁰(Ω¹_P2)=0なので、

`H⁰(Ω¹_X(2F))=H⁰(O_X)=C`。

従って各fixed projectionについて、(1)が与えるπ_Xのprojective sourceは高々一点。
C_S divisorはνから一意に定まり、自由な追加parameterとして数えていない。
ambient πもπ_Xから高々一つに定まる。
実際、二つのambient bivectorがDのdense smooth locusで等しければ差はE(−5)のglobal sectionだが、
`H⁰(P4,E(−5))=H⁰(P4,Ω²)=0` である。

## 6. 全projection parameterと実stabilizerによるdimension bound

完全H⁰(X,L)は八次元で、projectionは五次元subspaceVとそのtarget projective frameで指定される。
frameのdimensionは24。
plane Pの二つのlinear equationsは、pullbackするとC_Pに沿って零となる二つのV-sectionを与える。
nodeでは各branchに一次、cuspではt²,t³なので唯一のbranchに少なくとも二次で零となる。
従ってC_P class2Fの両caseを通じて

`H⁰(L(−C_P))⊂V`

となる。
実際、両scrollで `h⁰(L−2F)=h⁰(E(−2))=2` であり、二つのplane equationは独立なので
この全二次元spaceを尽くす。

C_P∈|2F|のchoiceはdimension二。
fixed C_Pに対し、それを含むVはGr(3,6)でdimension九。
finite birational、正しい導手型、πがambientへdescendする等の条件はこのsupersourceを狭めるだけ。
従ってprojectionの全parameter upper boundは `2+9=11`。
projective π_X choiceは前節でdimension零以下である。

actual polarized automorphism subgroupは
bundle automorphism modulo scalarとbase PGL2からなる。
X_113では `h⁰(EndE)=11` なのでdimension `11−1+3=13`、
X_122では `h⁰(EndE)=9` なのでdimension `9−1+3=11`。
これらは抽象的parameterの重複ではなく、X上の実automorphismである。
target frameも含むparameter source上でこの作用はfreeである。
もしautomorphismがV及びそのtarget frameをprojectively固定すれば、νと合成して同じfinite birational mapを与え、
Xのdense open上でidentity、従ってX全体でidentityとなる。
よって各image fiberはこのdimension13又は11のorbitを含む。

全PGL5像のdimensionはそれぞれ

`24+11−13=22`, `24+11−11=24`

以下となる。
repeated-rootのC_P、rank零又はrank二plane、special projection centresも全て同じparameter upper boundに含まれる。
有限unionなので未知のdimension≥25 componentは生じない。命題が従う。

## 7. 残るdegree五分岐

length formulaでdegree五rank零Sはδ_S=1である。
本noteは残余planeを持つg=0の全場合を扱った。
残るg=1では二次元の非正規支持はSだけとなり、Dの正規化のsectional genusも一。
その全normal polarized threefold分類とPoisson導手の橋は本noteでは証明していない。
全rank零degree五又は全Gの完了とはしない。
