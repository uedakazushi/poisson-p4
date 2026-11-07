# smooth cubic scroll：負のdirectrixを零因子が含む必要条件の独立検証

2026-11-07、別作者 `/root/quartic_review`。**PASS（純証明）**。
producer `CUBIC_SCROLL_NORMAL_ROOT_DIRECTRIX_REDUCTION.md`、SHA256
`331c45892887bc2bc5fb80f201252f71fea42bfaec056e1560e43bc800f714f2` を全文読解した。
新しいsource/minor/Jacobi/GB/jet計算はない。既21源の数値dimension結果をこの命題の証明に使わない。

対象はsmooth scroll S=F1=S(1,2)⊂P4、H=C0+2f、actualD被約既約quintic、
S⊂SingD、restriction πSのgeneric rank2。
R=divπS∈|2C0+3f| の**全可約・非被約**を含む。
結論は `C0⊂R` の必要条件だけであり、全G上界・全空性・E membershipではない。

## 全域normal bundleとprimitive下界

三quadric q0,q1,q2と二Hilbert–Burch syzygyを直接展開すると、
scroll parametrization `(su,tu,s²v,stv,t²v)` 上でsyzygyのimageは
`(tv,−sv,u)` が生成する一つのline subbundleになる。
成分はclass L=C0+fのsectionで、同時に零とならない。
従って全域で

```
0→O(-L)→O³→N*(2H)→0,
N*(2H)=β*Q, β=[tv:−sv:u]:F1→P2.
```

βはC0をp=[0:0:1]へ縮約するblowdown、β*O(1)=L。
単なるc1/c2一致からbundle同型を推測しない。
ruling fへの制限はN|f=O1⊕O2、directrixはN|C0=O2²。
後者は注入 `(s²,st,t²):O(-1)→O1³` の二linear syzygyからも確認できる。

primitive factor ℓ:N→A、A=aC0+bfでは、一般rulingとC0へのnonzero restrictionから
`a≥1,b−a≥2`。C0上全写像零ならC0というcommon divisorial zeroなのでprimitiveでない。
一方

```
Hom(N,A)=H0(P2,Q(b−4)⊗I_p^(b−a−2))
```

はpullback/twistとβ_*O(-mC0)=I_p^mのprojection formulaで全域に成立。
Euler sequenceとH1(OP2(d))=0からH0(Q(d))=0 for d<0、従ってb≥4。
codimension2 basepointはこの証明から除去していない。

constant-vector pのQ-sectionはpでsimple zeroだけを持つ。
blowup各chartでexceptional factorを割ると(1,t)/(s,1)、従ってnowhere-zero sectionを得て

```
0→O(2C0+3f)→N→A0=O(C0+4f)→0.
```

a1なら左kernelからAへのmapはruling degree−1なので零。
全mapはこのquotientとH0(O((b−4)f))の積で、b>4はfiber divisorial zeroを持つ。
よって**a1 primitive rootはA0、ℓ0のconstant倍だけ**。

## ambient dq根の全排除

A=2Hの全HomはH0(Q)=C³、ambient quadric netの全normal derivativeそのもの。
c0又はc1非零のQcはrank4、c0=c1=0のq2はC0上でdq2が零なので非primitive。
後者のsaturationはちょうどℓ0であり、primitive A=2H例として残さない。
rank4性は四変数minorに加えてkernel vector
`(-c1c2,c0c2,c1²,-c0c1,c0²)` でも確認できる。

rank4 Qがnormal quadraticの一根ならgeneric smooth S⊂Qでnormal座標(q,z)を取り、
D∈(q,z)²のquadratic項にqが因子として入るためD|Qは3Sを含む。
rank1 normal quadraticのdouble rootにも同じ含意が適用される。
Qのvertexは点なのでgeneric Weil multiplicityを変えない。
Cl(Q)=Z²、H=(1,1)、integral cubic scrollは(1,2)又は(2,1)。
effective classは二ruling一般直線の交点数から両係数非負。
(0,3)/(3,0)は一ruling familyのplanesの和でintegral nonplane scrollにならない。
D|Q≠0のresidual classは(2,−1)/(−1,2)で非effective。
よってQ|Dが必要となるがD既約に反する。
**A=2Hのprimitive根は両normal ranksで全排除**。

## 全 normal rank2：八sumclassと五pair

nonplane scrollのsecant argumentはnormal multiplicity3以上をdegree5と矛盾させ、一般位数は2。
Weinstein generic rank2分解とD∈I²からconormal bracketはabelian。
G/πSはflat parallel normal quadratic、従って
`discG=c(πS)²`、G≠0。
ここはRのsmoothness・被約性を使わない。

c≠0はC(S)上でsplitする。primitive global factorsをℓi:N→Aiと取り

```
G=γℓ1ℓ2, V=divγ, Δ=div(ℓ1∧ℓ2),
V+Δ=R, V=5H−A1−A2, Δ=A1+A2−detN.
```

height1でprimitiveなのでγにdivisorial poleはなく、smooth SのHartogsでcodim2も延長できる。
Ai=aiC0+bif、a=a1+a2,b=b1+b2に対しeffective coneと下界は
`3≤a≤5,8≤b≤10,b−a≥4`。
全八sumは
`(3,8),(3,9),(3,10),(4,8),(4,9),(4,10),(5,9),(5,10)`。
V/ΔがC0を含まないと仮定する誤ったnef条件を加えない。

V0はRが非被約でも不可能。
E=P(N)上のclosed residueαのdivisorial polesは二horizontal prime Γ1,Γ2だけ。
**Δでrootが衝突してもそれはΓ1∩Γ2というcodim2集合で、余分なvertical pole divisorではない**。
γにzeroが無いのでstrictD∩Eはgeneric vertical divisorを含まず、αもそこでregular。
二Γ上のresiduesはnonzero λ,−λ。
一般ample complete-intersection curvesは全codim2 bad lociを避け、曲線留数和から
`λ(Γ1−Γ2)≡0`、従ってAiが数値同値。
Pic(F1)=Z²だからAi同一、しかしA1+A2=5HはC0係数五がoddで矛盾。
全点SNC又はordinary-root collisionの不在を仮定していない。

このV0排除、a1唯一根、ambient根排除を全八sumへ適用すると、残るunordered pairsは正確に

| primitive pair | V | Δ |
|---|---|---|
| (1,4)+(2,5) | 2C0+f | 2f |
| (1,4)+(2,6) | 2C0 | 3f |
| (1,4)+(3,5) | C0+f | C0+2f |
| (1,4)+(3,6) | C0 | C0+3f |
| (2,5)+(2,5) | C0 | C0+3f |

第三行以外はV·C0<0又はVの明示形からC0を固定成分として持つ。
第三行V=Lがprimeならγのvertical primeはreducedで**generic poleはsimple**。
R/Δが同じprimeを更に含む可能性も保持するが、一般fiberのlinear rootsは非零なので
vertical pole orderはγの1であり、higher poleを仮定していない。
closednessと一般曲線留数和は
`λ(A1−A2)+μV=0` in Pic(F1)⊗Cを与える。
Ai差=−2C0−f、V=C0+fは独立、λ≠0なので矛盾。
従ってV非prime、|C0+f|の非primeはC0+f0だけであり、やはりC0⊂V⊂R。

## 全 normal rank1と結論

c0ならG=γℓ²、V=5H−2A effectiveかつRに支持。
下界から1≤a≤2,4≤b≤5。
(1,5)は非primitive、(2,4)はambient根で排除。
残A=(1,4),(2,5)はV=3C0+2f,C0で、両方ともC0固定因子。
G0は一般位数2と矛盾するので別の隠れた補集合ではない。
これでnormal rank1/2と全Rの可約・非被約を尽くし、**C0⊂R** が従う。

R|C0=O_P1(1)なので、restrictionπSをC0に制限したsectionの二linear係数が零という
全域必要条件を元21源に戻してよい。
その二係数を加える次のfinitejobは支持するが、今回その数値rank/G上界を確認したとはしない。
smooth scroll以外、twisted cubic cone、全degree3 rank2支持、全G/Eの完了は未解決として保持する。
