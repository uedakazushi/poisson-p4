# Canonical index二正規化の三重曲線代数：復元した部分還元

2026-11-07、G、producer plane_collision_review。
これはVM交代前の未pushノートの**数学内容を会話の引継ぎから復元した別identity**である。
旧本文のbyte単位復元ではない。旧expected SHA256は
`54474918c742f86be8e6378580bb78e05c66dcc40f1b583ae411cd67e3b4e2f6`。
旧manifest expected SHA256は`65ade9bef79550d548a7101bd94344bee15c4694c7e0134f846434f80e1f00c5`。
新CAS無し、独立review前の部分証明で、これだけで全G上界又は全E被覆を主張しない。
後続のwhole証明は本復元ノートを依存にしない。

## 1. 入力と全scheme

入力は `RANK_TWO_QUINTIC_CANONICAL_INDEX_TWO_REDUCTION.md`。
Xはnormal canonical Gorenstein、−KX=2L、L ample/GG、L³=5、
ν:X→Dは既約被約quinticの有限正規化、Sはintegral非退化degree五の導手支持。
rank二nodeの二prime又はcuspの一primeに対しΓ∼Lである。
全normalization moduleは

```
M=[[C,A,B],[A,u0,u1],[B,u1,u2]], q=u0*u2-u1²,
F=ν_*OX=coker M,  F0=O⊕O(-1)^2,
I_S=(q, B*u1-A*u2, A*u1-B*u0).
```

u0,u1,u2は独立、V=(u0=u1=u2=0)はP1、全M entriesはI_Vに属する。
T=ν^{-1}Vはschemeとしてfinite flat三重被覆で
`ν_*OT=OV⊕OV(-1)^2`、degree三、arithmetically CM、pa零、P3をspanする。
Tのreduced性、Γのnormal性、Xのterminal性は仮定しない。

## 2. SとΓのHilbert式、およびV上の差

vanishingと−KX=2Lから
`Hilb(R_X)=(1+3t+t²)/(1-t)^4`。
conductor moduleはF(-2)であるから、
`Hilb(R_D)=(1-t^5)/(1-t)^5`を使い
`Hilb(R_S)=(1+2t+2t²)/(1-t)^3`を得る。
SはACM、一般hyperplane曲線の算術種数は二。

I_S⊂I_V²である。Vのgeneric点におけるSのtransverse一次元CM局所環Rは
embedding dimension三、従ってmultiplicity e≥3。
normalization Bとのδは
`length(B/R)≥dim_K(B/m_R B)-1=e-1≥2`。
一般hyperplane曲線の算術種数二より、Vに沿うδは正確に二、正規化曲線は種数零で、
他の非正規曲線は無い。SはS2なのでS\Vはnormalである。

各Γ∈|L|について
`Hilb(R_Γ)=(1+3t+t²)/(1-t)^3`。
finite birational Γ→SのcokernelはHilbert式`t/(1-t)^2`、CM pure curveでVに支持される。
generic length一なのでI_Vで殺され、P1上のrank一CM sheafとして正確にO_V(-1)である。
従って各ΓもVの上で非正規であり、その一般hyperplane曲線のδは一。

## 3. Generic三次代数上のmark

K=C(V)、A=(ν_*OT)_Kとする。dim_K A=3。
Γのdefining section sのTへのrestrictionをaとする。
Γ→Sのcokernel O_V(-1)をKへ制限すると、A/(a)からKへの非零線形写像がありunitを殺す。
unitそのものはA/(a)で非零であるから、dim A/(a)≥2、multiplication by aのrankは≤1。
a=0ならs∈span(u0,u1,u2)⊂ν^*H0(P4,O1)となり、Γの像Sがhyperplaneに入る。
非退化性からこれは不可能なのでrankは正確に一である。
nodeでは二mark a,bが非零rank一でab=0、cuspではaが非零rank一でa²=0。

## 4. Reducedの場合とnonreducedの場合

T reducedならAはreduced separable cubic K-algebra。
cubic fieldでは非零multiplication rank三、quadratic field×Kではrank一元は唯一のK因子にしか支持できない。
二つの非零rank一元が積零となるnodeではA=K³だけが残る。
Tの三つの成分はそれぞれVへdegree一、finite birationalに写るので各P1である。
TはCM degree三、pa零でP3をspanするから、三直線のchain又は非共面concurrent tripleとなる。
conic+line型はこの二mark条件を満たさない。

T nonreducedならgeneric Aもnonreducedである。
そうでなければnilradicalはOV-torsionとなり、OV-flatなOTへの単射に反する。
三次元commutative nonreduced K-algebraの全型は

```
K × K[ε]/(ε²),  K[ε]/(ε³),  K[ε,η]/(ε,η)².
```

residue field degree≥2を持つ非自明nilpotent因子はそれだけでdimension≥4なので現れない。
cuspのrank一nilpotentは各々(0,ε)、ε²、非零nilvectorである。
このgeneric代数分類だけからglobal OV-algebraやXのextensionを分類したとは言わない。

## 5. 残余と証拠の区別

これは全三次代数のgeneric型とΓ→Sの差の部分還元である。
global flat familyのnormal form、Xのnormal方向extension、Poissonのambient descent、
source次元上界はこの時点で未完了だった。
新CASや全源計算を実施したという記録ではない。
後続whole証明で不要になっても、本内容の旧partial scopeを遡って変更しない。
