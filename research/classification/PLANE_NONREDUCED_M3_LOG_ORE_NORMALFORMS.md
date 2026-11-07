# 六残constant-root源の全module有限化：十whole lineと十二fixed

2026-11-07、主義務G。producer純証明。
前note `PLANE_NONREDUCED_M3_CONSTANT_ROOTS_ORE_GEOMETRY.md` の三閉鎖と
E1 EXACT_B/TRIPLE_LINEARの還元を使う。新module算術・全constant-root G閉鎖はまだ主張しない。
本noteでは残LOG e二の全supportを尽くし、有限moduleリストを確定する。

## 1. 記法と全Q-space

y=(a,b,z,w)、Π=Q+∂c∧W、全W=W0+W十八。
LOG e二では全α=tに対しW0=(a²,−2ab,−2az,3taw)。
Fは全normal quadratic。d=3t+2と置く。
d非零かつt≠−8/3では全Qは五cohom r,B,C,D,EとFで

```
Qab=λBbw+Caw+kDbz+F,
Qaz=λBzw+kDz², Qaw=νBw²+ℓDzw,
Qbz=−rz²+σCzw, Qbw=(E−r)zw+ρCw², Qzw=0,
λ=(4−3t)/d, ν=(2−6t)/d, k=−6/d, ℓ=(3t−4)/d,
σ=(1−3t)/d, ρ=(−6t−1)/d.
```

generic completeQ-space、元operatorの全壁specializationをreceiptで照合した。
この八dim familyはt=1/3,4/3でも全rank八を保つ。
Bのλ/νは同時に零にならず、Dのbz項、Cのaw項、rのz²項、Eのzw項から独立性が分かる。
t=−8/3の新一方向、d零は以下で別に復活する。

strict Jacobiの全reduced supportは正確に

`rD=0, BC=0, dB(r−E)+(d+3)CD=0`

である。係数にあるt,3t−1又は3t−4を個々に反転せず、各対のgcd一から全tで同じ三式を得る。
actual D3は `a²z²w*(dr−3E)`。従って対象m三ではdr−3E非零。

## 2. generic全support

B非零ならC零、E=r、r非零、D零となり、t≠1/3。
そのQをB型と呼ぶ。
B零・C非零では、t≠−5/3ならD零、Qはb非依存の∂b定数因子となる。
t=−5/3でD非零ならr零、E非零で、下記critical CD型へ入る。
B=C零ではD零なら同じ定数因子、D非零ならr零、E非零でD型となる。
従ってこれで全zero/nonzero補集合を尽くした。
定数因子の全actual Dは既essential三八成分へ入り、全像二十三以下。

D型は非零cohom二係数を実diagonalsで一にし、

```
Qab=bz+F, Qaz=z², Qaw=Rzw, Qbz=0, Qbw=zw,
R=(4−3t)/6
```

となる。central coboundaryはspan(zw)なのでF=εz²+ηw²、ε,η∈{0,1}。
全tの四whole lineを保持する。

B型ではνBが非零なのでそれと−rを実diagonalsで一へ戻し、

```
Qab=Λbw+F, Qaz=Λzw, Qaw=w², Qbz=z², Qbw=0,
Λ=(4−3t)/(2−6t)
```

となる。t≠1/3。central coboundaryはspan(z²,(1−Λ)w²)。
t≠−2/3ならF=εzw、ε∈{0,1}。二whole lineを保持する。
t=−2/3で消えるw² coboundaryは消去せず、下記scalar B四fixedへ残す。

## 3. scalar t=−2/3全source

元全kernelから直接、五cohom A,r,B,C,EとFの全Qは

```
Qab=−Abz−Bbw+F,
Qaz=−Az²−Bzw, Qaw=−Azw−Bw²,
Qbz=−rz²−Czw, Qbw=−rzw−Cw²+Ezw, Qzw=0.
```

全strict Jacobiは `Ar=0, BC=0, AC+Br−BE=0`、D3非零条件はE非零。
A=B零なら定数因子。
A非零ならr=B=C零となりD型のt=−2/3へ入る。
B非零ならA=C零、r=Eでscalar B型へ入る。
A,B両非零はE非零に矛盾する。

scalar BはQab=bw+F,Qaz=zw,Qaw=w²,Qbz=z²,Qbw=Qzw零。
central coboundaryはspan(z²)だけで、F=εzw+ηw²、ε,η∈{0,1}の四fixedを保持する。
全四でW0=(a²,−2ab,−2az,−2aw)。

## 4. 二つのcritical支持

t=−5/3のCD非零枝はr=B零、C,D,E非零。
実diagonalsで三係数を一にして

```
Qab=aw+2bz+F, Qaz=2z², Qaw=3zw,
Qbz=−2zw, Qbw=−3w²+zw, Qzw=0.
```

central coboundaryはspan(zw,w²)、従ってF=εz²の二fixed。
W0=(a²,−2ab,−2az,−5aw)。

t=−8/3の全九dim Q-spaceには新係数Aがあり、元receiptの座標をr,A,B,C,D,Eとすると

```
Qab=−2Bbw+Caw+Dbz+F,
Qaz=−2Bzw+Dz², Qaw=Az²−3Bw²+2Dzw,
Qbz=−rz²−3Czw/2, Qbw=(E−r)zw−5Cw²/2.
```

全strict Jacobiは
`rD−3AC/2=0, 2B(E−r)−CD=0, BC=0, A(r+E)=0`。
D3非零は2r+E非零。
A零は§2のgeneric eight-spaceへ戻る。
A非零ならE=−r、r非零。三式から
D=3AC/(2r)、B=−3AC²/(8r²)、さらにBC=0よりC零、従ってB=D零。
実diagonalsでA=r一、

```
Qab=F, Qaz=0, Qaw=z², Qbz=−z², Qbw=−2zw, Qzw=0.
```

central coboundaryはspan(z²,zw)、F=εw²の二fixedを保持。
W0=(a²,−2ab,−2az,−8aw)。
これでextra Aの零/非零も尽くした。

## 5. 実normalizer、central全source、逆

LOG first symbolを保つ全対角変換の十分な部分群は
`(a',b',c',z',w')=(sa*a,sb*b,sc*c,sz*z,sw*w)` とtensor scale sa/sc。
五scaleは非零、a/b比とnormal二scaleは自由である。
従って各branchの指定した二又は三非零cohomを同時に一へ戻せる。
EXACT_B e一とTRIPLE_LINEARではsb=saという追加条件があるが、
先行noteで使う一非零cohomと二central零/非零の正規化には十分な自由度が残る。

actual a→a+u(z,w), b→b+v(z,w)でu,v任意linearとすると、
Qabのcentral Fの変化は

`δF=Xa(v)−Xb(u)−(∂a Qab)u−(∂b Qab)v`。

ここでXa=Qaz∂z+Qaw∂w、Xb=Qbz∂z+Qbw∂w。
各normal formへ代入すると上記central imageを正確に得る。
単にFを任意の値へ固定しておらず、image外の全coefficientsを残した。
残diagonal scalesは各residual coefficientを独立に零又は一へ直す。
例えばD型では固定cohomを保つ条件sa=sb, sz=(sa/sc)saの下、
Fのz²/w²係数は1/σとσsa²/sw²倍、σ=sa/scなので独立。
scalar Bではsw=σsa,sz=σsbの下、zw/w²係数は1/σとsb/(σsa)倍で独立。
他一coefficient残余も自由σ又はnormal比で正規化出来る。

全inverseはminus shearとnonzero diagonal inverse。
base-normal shearはP0のweight零を保ち、変化したWはW0+全十八sourceへ戻る。
cは移さないのでc非依存性も保持する。
残Euler gauge W→W+ℓE_yはmodule/actual projective tensorを保つ。
従って全higher-source包含の逆橋がある。

## 6. 一回の全module義務と上界

全未閉sourceは、定数因子又は既factor空枝を除き、次の有限全支持リストへ入る。

- E1 EXACT_B：εz²+ηw²、四whole t line。
- LOG D型：εz²+ηw²、四whole t line。
- LOG B型：εzw、二whole t line。t=1/3はleading零、−2/3は次fixedへ。
- scalar B：εzw+ηw²、四fixed。
- critical CD：εz²、二fixed。
- critical A：εw²、二fixed。
- TRIPLE_LINEAR：εz²+ηw²、四fixed。

合計十whole lineと十二fixed、全W十八を保持する。
各strict QはPoissonなのでfull350 projective JacobiはWにaffine linear。
一回のbounded finite batchで全Wmodule/augmentedrank/particular/completekernel、
actualD全affine係数のproper factor、又は弱いdimensionを確認する。
全generic minorとfactor-specialization分母の全根を元matrixへ戻す。

dimension fallbackでは固定tのQ-source dimension upperを、E1は4、LOG D/Bは5、
scalar Bは5、CDは6、critical Aは5、TRIPLE_LINEARは5と取れる。
central零条件一つにつき一減出来る。
固定plane first-symbol source upper BはE1で9、LOG whole lineで8、
LOG fixedtで7、TRIPLE_LINEARで7。
全projective像は `B+dim(Q-source)+dim(Wkernel)+5` 以下。
全actual normalizationがsourceを全moduleへ送るので、canonical点のWkernel次元を全orbitへ戻せる。
fixedt wallではBを一減し、上界一般値を特殊値へ外挿しない。
module不整合/properDfactorならその全枝は空。

全二十二moduleの結論が未確定の間はconstant-root全Gを閉じない。
他primitive degrees/generic triple/m二/reducedRは別残余のままである。
