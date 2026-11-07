# E6 quarticとcuspを通るtransverse line：全族はrank一cofactorへ

2026-11-07、PRIMARY E。producer手証明、独立確認前。
H+既約quarticの九profile入口のうち、FがE6 cuspを持ち、Hがその点を通って
cuspidal tangentと異なるlineである全marked branchを扱う。
結論は、二次Saito frameを持つ全点のcanonical H-cofactorがrank一になること。
従ってtype6 nil22-imageの残余にはこのbranchは存在しない。
E6 cusp tangentのr1配置、cuspを通らないH、他のprofileは含めない。

## 1. 全marked quarticとτ12

H=h、cusp p=[0:0:1]、そのtangentをy=0に取る。
E6のmultiplicity三を用いると全quarticは

```
F=y³w+a h⁴+b h³y+c h²y²+d h y³+e y⁴.
```

a!=0である。a=0ならyがFを割り、integralityにも反する。
wの可逆shearでd,eを消せる。H|Fはy³wなので、cuspでI=3、他の一点でI=1。
Fはdelta三を使い切り他の特異点を持たない。

w=1でg=hFを取り、weights(h,y)=(3,4)とする。

```
g0=h*y³+a*h⁵, wt(g0)=15,
g=g0+b*h⁴*y+c*h³*y²,
wt(h⁴*y)=16, wt(h³*y²)=17.
```

二partialのinitial formsはregular sequence、weights十二と十一。
weighted Milnor algebraのlengthは十一、topweight十六である。
weight十六部分は(y⁴,h⁴y)をbasisとし、Jacobian idealの同weight部分は
`y*(g0)_h=y⁴+5a h⁴y`だけなので、h⁴yのclassは非零。
initial formsがregular sequenceであるため、gのfiltered Milnor algebraもこのassociated gradedを持つ。
weightedEuler identityにより

```
-15[g]=b*[h⁴y]+2c*[h³y²].
```

後者はtopweightより高く零。従ってcusp側tau十一はb=0と同値。
もう一つのH交点はsmooth transverse nodeでtau一なので、全tau(D)=12は正確にb=0。
weighted leadingだけでquasihomogeneityを推測せず、socle coefficientを全familyで比較した。

## 2. cの零補集合を含む全正規化

F=y³w+a h4+c h2y2、a!=0を考える。
対角変換でa=1に直し、残るcは零又は非零の二caseへ移る。
非零ならh,y,wの非零scaleでc=1に出来る。
元のa,cへは逆diagonalで戻れ、未知coefficientの境界c=0は次に別に扱う。

## 3. 二次logarithmic anchorを全係数で解く

h=1でF=y³w+1+c y²。全quadratic field(U1,U2)と全affine linear ellに
`U1*F_y+U2*F_w=ell*F`を課す。
両U成分各六係数、ell全三係数を比較し、選んだ二導分だけに制限しない。

c=0では全解が

```
(U1,U2)=(u+a*y+b*w)*(y,-3w)
```

となる。全三係数を持つが、任意の二本のdeterminantは零。
従って二次Saito basisを持たず、元のfull-rank cotangent sourceの仮定を満たさない。

c=1では全解は次の二本の定数係数線形包。

```
U=(y²,-2-3yw),
V=(yw-2/3,-4/3-(4/3)yw-3w²).
```

直接

```
det(U,V)=-(4/3)F,
U(F)=0,
V(F)=-(4/3)y*F,
[U,V]=-5w*U+y*V
```

である。従ってhomogenized E,U,Vは真の全quadratic Saito frameとなる。
全六磁気Rに対し

```
P=partial_q1 wedge U+partial_q2 wedge V
  +(-5w*q1+y*q2+R)*partial_q1 wedge partial_q2
```

が実Poisson sourceである。元の全sourceframeは定数GL2 fiber変換で、全magneticは
可逆な二次形式pullbackで戻せる。

## 4. 全六磁気でのcanonical cofactor

一般のcotangent sourceで[U,V]=alpha U+beta V、D=det(U,V)とすると

```
alpha=div V-V(logD), beta=U(logD)-div U,
div(P)= -U(logD)*partial_q1-V(logD)*partial_q2.
```

任意のRはbaseだけに依存するので、この公式に磁気係数は現れない。
上の二本ではdiv U=−y、div V=−(4/3)y−5wであり、

```
div(P)=(4/3)y*partial_q2.
```

従って全sourceでcanonical cofactorは平方零・rank一、affine constant部零である。
rank二nil22の入口へは入らない。
全Hquarticとして見ればnil211・translation零の既済whole sourceに接続できるが、
そのMINIMAL_M帰属は既存の全源定理への依存であり、今回新しいambient点を計算していない。

本結果は全marked coefficients、c零補集合、全六magneticを保持したbranch排除である。
全Hquartic・全Eの完成や他profileへの自動拡張は主張しない。
