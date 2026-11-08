# 二つの weak wall の実二十三軌道と既二十三親の排除

2026-11-08、PRIMARY E、producer 純証明・独立 review 前。
新しい tangent/Jacobi/operator/rank/minor 計算は無い。
入力は `PLANE_M3_TRIPLE_WEAK_WALLS_TWO_ACTUAL_ORBITS.md` の
full W六 transitivity／Dmap injectivity と二つの particular tensor。
ここでは actual divisor stabilizer の elementary coefficient 比較だけを使う。

## 1. divisor stabilizer の全可能 Lie field

二つの D を共通に

`D=K ab z³+L a w³z+M b w²z²−K c z⁴+E w⁵`

と書く。A型(τ4/21)は L零、B型(τ1/12)は L非零、他 K,M,E は全て非零。
既証明の reduced singular support は唯一の plane S={z=w=0}。
generic normal leading は `Kabz³`。
従って connected projective stabilizer は S、normal triple root z、
S上の unordered pair of lines a零/b零を保つ。
連結群はこの二 lines を個別に保つので、可能な全 linear Lie field は

```
Xa=αa+A z+B w,
Xb=βb+C z+D w,
Xc=γc+E1 a+F1 b+G z+H w,
Xz=λz, Xw=μw+νz.
```

normal rows に base 項を落としたのは S保存、a/bにcを混ぜないのは
base UFD identity a'b'∝ab からであり、選んだ座標部分群だけを数えていない。

## 2. off-diagonal coefficients の全比較

X(D)=κD を比較する。
B型では monomial a w²z² の coefficient が3Lνなのでν零。
A型では w⁴z の coefficient が5Eνなので同じくν零。
次に b w z³ は KB、a w z³ は KD なので B=D零。
w²z³ は MC なので C零。

B型では w³z² の係数LAが残るためA零。
その後az⁴,bz⁴,z⁵,wz⁴から E1=F1=G=H零。
従って B型に off-diagonal connected stabilizer は無い。

A型ではAだけが残り、bz⁴の係数からF1=A、
残monomialsからE1=G=H零。
唯一の unipotent direction は

`N_D=z∂a+b∂c`。

これは literal D を正確に保つ。

両型の diagonal 条件は、二自由 weightsを r=weight(z),q=weight(w) として

`weight(a)=2q−r, weight(b)=3q−2r, weight(c)=5q−4r`。

共通 Euler scale を除いた projective torus は一次元。
例えば weights(c,a,b,z,w)=(0,3,2,5,4)。
従って connected divisor stabilizer は A型で二次元、B型で一次元。

## 3. tensor stabilizer は両方とも一次元

この torus は normalized tensor scale を入れて literal Qτ と元 P0 を保つ。
weighted D は degree二十なので actual D0 も保つ。
群作用後の W を canonical Euler gauge に戻すと同じ full W六 source、同じ actual D0 となる。
その map の injectivity により、particular projective tensor は torusで固定される。

A型の unipotent N_D は tensor を固定しない。
実際 `a'=a+t z,c'=c+t b` の下で Qのcoreは

`Qab→Qab+(1−3τ)t z²`、他core entriesは不変。

τ4/21では1−3τ非零。
この difference は Euler redundancy ではない。
なぜなら `E_y∧(z²∂a∧∂b)` は非零であり、Euler bivectorならこの wedge は零だからである。
unipotent がprojective tensorを scalar倍する可能性もない：
unipotentのcharacterは自明、又は literal D不変と tensor-square のD変換からそのscalarは一である。
既に tensorを保つ torus方向を足しても、この nonzero unipotent quotient は消えない。

従って両 tensor の connected projective stabilizer は一次元、
**両 actual PGL5 orbit の dimension は24−1=23**。
char零のalgebraic stabilizerはsmoothなので、hidden infinitesimal groupを数え漏らしていない。

## 4. 二つの既二十三軌道 closure には入らない

known K22/nil22-long canonical orbit は二十三次元で、その generic actual divisor は
三つの distinct singular planesを持つ。
K4a canonical orbitも既ambient証明書で二十三次元、actual divisorは
{u=v=0} と {u=a=0} の二 distinct planesを持つ。
本二 weak orbit の actual divisor は一つの planeのみなので、どちらのorbitとも同じではない。

irreducible orbit closure の次元が二十三で、その dense orbit が二十三次元なら、
異なる二十三次元orbitを内部に持てない。
仮に持てば両方がそのirreduciblevarietyのdense openとなり、互いに交わる必要があるからである。
したがって本 weak points は **K22の既二十三成分又はK4aの二十三 orbit closure に含まれない**。

この議論は高次元 parent を排除しない。
特に K4a は smooth unique二十四親の境界と独立に証明されたため、
その二十四親、MPS二十四／二十五 component、又は未知 component への包含は未解決。
二つの orbit が full component であるとも、成分数を増やすとも主張しない。
