# μ5/r1：A5及びA4+A1全quartic familyのmarked lineを三orbitへ

2026-11-07、PRIMARY E。producer手証明、独立確認前。
tau(HF)=12の数値入口で、mu(F)=5ならr1、かつD=HFの全singularityはquasihomogeneousである。
ここではA5及びA4+A1の全global quartic familyを正規化し、marked Hの全可能性を列挙する。
二次Saito source/cofactor・全磁気・component membershipは別の一括検証に残す。

## 1. A5：一parameterを保持してからhyperflexを解く

正規化P1上、A5の二branch逆像をs=0,t=0に取る。
共通tangent hのpullbackはs²t²。transverse line yは各逆像で一次なので、
P1とbaseの可逆scale及びh-shearでy=st(s²+t²)へ直せる。
第三sectionのh,y成分を除けば

```
w=A s4+B s³t+C st³+D t4.
```

二branchのh/y²の二次leading coefficientはそれぞれD,A。
接触が少なくとも三であるためA=D!=0。
wのscaleとy-shearにより全quarticは

```
[h:y:w]=[s²t²:st(s²+t²):s4+t4+k(s³t-st³)]
```

へ移る。三次branch接触の差は2kなのでA5であることはk!=0と同値。
k=0ではmapはs↔tで二対一となり、quartic birational normalizationではない。
従って一parameter kをこの段階で捨てていない。

HがA5点を通ってintersection四を持つなら、Hは共通tangent hである。
しかしDの三smooth branchのpairwise contactは(3,2,2)。
weighted homogeneous三smooth branchにはこの接触表は無いのでtau12条件に反する。
従ってHはsmooth点でのhyperflexである。

general lineのaffine pullbackは

```
C*r4+(Ck+B)*r3+A*r2+(B-Ck)*r+C.
```

A5点を避けるのでC!=0、全intersection四が一smooth点ならC(r-u)^4となる。
constantからu4=1、r3/r係数からk=2(u³-u)。
u=±1はk0で除外、u=±iはk=∓4iでありbranch交換で同一marked orbit。
r→i*rと可逆base scaleにより、唯一のrational representativeは

```
F=(hw-y²-2h²)²-16h²y²-64h4,
H=w+6h.
```

対応するnormalizationは
`[h:y:w]=[s²t²:st(s²-t²):s4+t4-4st(s²+t²)]`。
Hのpullbackは(s-t)^4、F|H=y4である。
Fの唯一A5とHのsmooth hyperflexがA7を作り、total tau=5+7=12。
Fはwに関してprimitiveで、discriminantの非平方因子y²+4h²により絶対既約。

## 2. A4+A1：全oneparameter family

A4 cusp逆像をs=0、nodeの二逆像をs/t=±1に取る。
A4 tangent sectionはh=s4、cuspとnodeを結ぶlineからy=s²t²を取れる。
第三sectionのh,y成分を除けばw=a s³t+b st³+c t4。
nodeの二逆像のimageが同じことからa+b=0。
A4条件とbase-point-free条件からa,c!=0。
従って全familyは

```
[h:y:w]=[s4:s²t²:k(s³t-st³)+t4], k!=0,
F_k=(hw-y²)²-k²hy(h-y)².
```

nodeは[1:1:1]であり、二tangentのdeterminantはkの非零倍なので全k!=0でordinary node。
cuspはA4で、delta三を使い切るので他のsingularityは無い。
Fはwのprimitive二次式、discriminantの非平方因子hyにより絶対既約である。

### A4 tangentはtau12に入らない

r1でHがA4点にあれば接線で、contactは四。
Dがquasihomogeneousなら、A4 branchとsmooth Hは同じweighted coordinatesでhomogeneousとなる。
A4のweights(2,5)において別のsmooth weighted branchとのintersectionは二又は五だけ。
四は不可能なので、このmarked caseはtau12条件から除外される。
node又はsmooth点をHの唯一交点とするcaseは次に保持する。

### nodeを通るH：一orbit

nodeでのcontact distributionは(3,1)又は(1,3)でなければならない。
general pullback

```
A*r4+kC*r3+B*r2-kC*r+C
```

を(r-1)^3(r+1)又は(r-1)(r+1)^3と比較するとk=±2だけ。
r→-rはkの符号を交換するので一marked orbit、代表は

```
F=(hw-y²)²-4hy(h-y)²,
H=h-w.
```

F|H=(h-y)^4。nodeに三本のsmooth branchでcontact(3,1,1)のD8が生じ、
他のA4と合わせてtotal tau=8+4=12である。

### smooth hyperflex H：一orbit

Hがnode,cuspを避けるならpullbackはC(r-u)^4、u!=0,±1。
constant/r3/r係数よりu4=1,u²=-1,k=∓4i。
二caseは同じmarked orbit。
r→i*rでrational representativeは

```
F=(hw-y²)²-16hy(h+y)²,
H=h+6y+w.
```

F|H=(h-y)^4で、intersection点はsmooth。
FのA4+A1及びHとのA7からtotal tau=4+1+7=12。

## 3. 全complementsと次の検証

A5共通tangent、A4 tangentの除外はtau12のquasihomogeneous条件による。
他のlineはr1のwhole根条件で尽くし、kの非零familyを一代表に先に固定していない。
以上でA5一orbit、A4+A1二orbit、計三つのwhole marked候補を得た。
tau12だけから全quadraticSaito frameを宣言せず、全18係数kernel・canonicalcofactor・全六magを
同じ有限source batchで検証する候補とする。
D5のsmooth hyperflex、A3+A2のsmooth hyperflex、2A2+A1の全marked lineは別残余として保持する。
