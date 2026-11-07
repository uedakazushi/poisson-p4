# degree六rank二支持：non-conical quintic退化因子では空

2026-11-07、主義務G。新producer幾何証明。
**命題。** rank四Poisson tensorの実Dが既約被約・non-conical quinticであるとする。
その非正規支持にdegree六integral surface Sがあり、π|Sの一般rankが二なら矛盾する。
Sのnormality・smoothness、全域ordinary-node性又はvertical divisor不存在を仮定しない。
conical Dのrank二枝は明示補集合として保持する。

## 1. 正規化の全有限型

`RANK_ZERO_SEXTIC_CONDUCTOR_EXCLUSION.md` §1–2のlength/Bertini/Höring論証は、
π|Sのrankを使用しない。
従ってδ_S=1、他に二次元非正規支持は無く、Dの正規化(X,L)は
degree五minimal-degree threefold in P7である。
その全型はS(1,1,3)、S(1,2,2)、S(0,1,4)、S(0,2,3)、S(0,0,5)。
最後の三つからのfinite linear projectionはconeを保つので、non-conical Dでは除かれる。
残るXはsmooth rank三scroll over P1である。

`PicX=ZL⊕ZF`, `K_X=−3L+3F`, `C_X∼3L−3F`。
δ一の横断一般型はnode又はcuspである。
cuspなら唯一のconductor prime ΓについてC_X=2Γとなり、PicのL係数のparityに反する。
従ってgeneric nodeだけが残り、C_Xはreducedである。

## 2. conductorが既約なら閉residue一形式に矛盾

nodeのnormalization branchはcompleted node ringの商なので、DのPoisson bracketは各branchにregularに誘導される。
branch idealはHamiltonian derivationで保存されるからπ_Xはconductorにtangentである。
したがってπ_Xは全codimension一でregularで、Ω_Xのpole divisorはC_Xである。
ambient inverse symplectic二形式のPoincaré residue

`α=contraction(π_X,Ω_X)`

はclosedで、C_Xに沿って高々simple poleを持つ。
rank零の場合のようにπ_Xがconductorで零であることは仮定しない。

C_Xが一つのprime Γなら、閉性によりgeneric log residueはconstant λ。
general complete-intersection curve Aをbad codimension二locusを避けて選ぶと、
α|Aの全residue和は `λ*(Γ·L²)=12λ`。
compact curve residue theoremからλ=0となる。
よってαはdivisorial poleを持たず、smooth X上でholomorphic一形式へ延びる。
Xはrational scrollなのでH⁰(Ω¹_X)=0。
しかし実Dが被約であることからπ_Dは一般rank二、αは非零である。矛盾。

## 3. 二つのconductor primeの全Pic classを除く

generic nodeは二つのbranchを持つので、conductorがreducibleなら
`C_X=Γ1+Γ2`、各Γ_iはSへgenerically degree一でfiniteに写る。
従って `L²Γ1=L²Γ2=6`。
各Γ_i=a_iL+b_iFと書く。
一般P2 fiberへrestrictするとa_i≥0で、a_i=0のprime effective divisorはfiber Fだけだが、
そのdegreeは一なのでここでは不可能。
従ってa_i≥1。
`a1+a2=3`、`5a_i+b_i=6` より、交換を除いて唯一の全可能性は

`Γ1∼L+F`, `Γ2∼2L−4F`。

しかし|2L−4F|にはprime divisorが無い。

X=S(1,1,3)でfiber座標をu,v,w（各degree1,1,3）とすると、
`H⁰(2L−4F)=H⁰(Sym²E(−4))` の全sectionは

`w*(a*u+b*v+c2(s,t)*w)`

である。全てw=0というquadric divisorを固定成分に持ち、integralでない。

X=S(1,2,2)では全sectionは、degree二blockの二座標v,wだけによるconstant binary quadric

`a*v²+b*v*w+c*w²`

である。複素数体上で二つのlinear factorへ分かれ、double factorの場合も含めprime divisorでない。
よって必要なΓ2は存在しない。

irreducible及びsplit conductorの両場合が除かれ、命題が従う。

## 4. 残る射程

本命題は全Sを含むがDのnon-conical性を仮定する。
conical Dにおいては正規化Xがsingular minimal-degree scrollである可能性を保持する。
その場合のrank二conical Sを、rank零専用の高次数cone排除で誤って除いていない。
degree五以下のrank二支持、plane restrictionのreducible/nonreduced cubic、及び全G/Eも別の残余である。
