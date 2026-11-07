# A4+A2 quarticの全marked line：r≤2を五配置へ尽くす

2026-11-07、PRIMARY E。producer手証明、独立確認前。
integral plane quartic FがA4及びA2だけを持つ場合を扱う。
全marked(H,F)のr=#H∩F≤2を、局所profileの一代表だけでなく全projective geometryで列挙する。
四つのr2配置の二次Saito基底・canonical cofactor・全磁気source帰属は、まだ計算確認していない。

## 1. 全quarticの実projective正規化

正規化はP1、二cusp逆像をs=0（A4）とt=0（A2）に取る。
A4の接線h=0は接触四である。
Puiseux(2,5)の接線contactは少なくとも四、Bezout上界四だから正確に四となる。
二cuspを結ぶline y=0は両cuspへtransverseでなければならない。
一方に接すれば交点総数が四を越えるためである。
従ってbase二sectionをscaleして

```
h=s⁴, y=s²t²
```

と取れる。第三sectionのh,y成分を可逆base shearで除けば
`w=a*s³t+b*st³+c*t⁴`。
t=0でA2 cuspとなるにはa=0、b!=0である。
s=0のcuspがbaseの点[0:0:1]に写るのでc!=0。
b,cの非零scaleとP1のs/t rescalingによって、全quarticは正確に

```
[h:y:w]=[s⁴:s²t²:12st³−6t⁴]
```

へ移る。逆は用いた非零scalarsとbase matrixの逆。
これは全marked quartic familyからの構成で、projective moduliを代表値へ勝手に固定していない。
実際にこの像は

```
F=(hw+6y²)²−144hy³
```

であり、既約・A4+A2であることは既済Ga9power3の入力でも確認される。

## 2. cuspを通る全line

一般line L=A h+B y+C wのpullbackは

```
A*s⁴+B*s²t²+12C*st³−6C*t⁴.
```

Fの正規化はcusps以外で単射なので、相異なるrootの数は相異なる交点の数に一致する。
cusp inverse自体も各一つなので、以下の全root比較でrを数えられる。

A4点を通るならC=0、pullbackはs²(A s²+B t²)。
A,Bがともに非零なら残る二rootは相異なり、r=3となる。
従ってr≤2はA=0のjoining line y=0、又はB=0のA4 tangent h=0だけ。

A2点を通るならA=0、pullbackはt²(B s²+12C st−6C t²)。
C=0は同じjoining line y=0。
B=0はA2 tangent w=0で、pullback6t³(2s−t)、contact(3,1)。
B,Cがともに非零なら二次因子がdouble rootとなる必要があり、判別式
`24C(6C+B)=0`からB=−6Cだけ。
line w−6y=0、そのpullbackは−6t²(s−t)²、contact(2,2)である。
このlineはA2へtransverseで、別のsmooth点に接する。

## 3. 両cuspを避ける全line

この場合A,C!=0。r=s/tと書くとpullbackは
`A*r⁴+B*r²+12C*r−6C`、零も∞もrootではない。
r≤2のroot partitionは(3,1)又は(2,2)、r1なら(4)である。

(3,1)でrootをu,vとするとr³係数零からv=−3u。

```
(r-u)³(r+3u)=r⁴−6u²r²+8u³r−3u⁴.
```

linear/constant係数を比較するとu=4/3。
唯一のlineは

```
h−(32/3)y+(128/81)w=0.
```

接触三のsmooth pointと、別のtransverse smooth pointを持つ。
(2,2)ならrootはu,−uで、polynomialはA(r²−u²)²。
linear coefficientは零となりC=0に反する。
(4)ならr³係数からu=0となり、cuspを避ける仮定に反する。
したがって他のsmooth bitangent/hyperflex lineは無い。

## 4. 全五配置と次の有限source課題

| marked H | r/contact | 配置 |
|---|---|---|
| h | 1/(4) | A4 tangent、既済Ga9power3 |
| y | 2/(2,2) | cusp joining line |
| w | 2/(3,1) | A2 tangent |
| w−6y | 2/(2,2) | A2 transverse＋smooth tangent |
| h−(32/3)y+(128/81)w | 2/(3,1) | smooth flex＋smooth transverse |

四つのr2配置の局所Dはそれぞれ

```
D7+D5,
A4+E7+A1,
A4+D5+A3,
A4+A2+A5+A1
```

であり、いずれもsimple、total tau=12である。
しかしtau12だけで全二次Saito frameや全Poisson sourceの存在を宣言しない。
前のE6transverseのc0のように、tau条件を満たしても二次anchorが全て比例する可能性がある。
各配置で全18quadratic gradient-syzygy係数と全六磁気を保持した一batchの検証を次に用意する。
cofactorがrank二nil22か、rank一既知sourceか、別のcanonical型かもその段階で確認する。
新component・既知closureへのmembershipを本noteだけで数えない。
