# Dependent H二全源：実shear、点ごとのP2一意性と全排除の独立監査

2026-11-08、PRIMARY E、別作者 `exceptional_coverage`。**純証明 PASS**。
producer `PLANE_DEPENDENT_H_SEMISIMPLE_WHOLE_SOURCE_EXCLUSION.md`（SHA5424d1fe…）と、
既独立六源receiptの全kernel・実shear・injective J2という入力を読解した。
新しいM1/J2、point、rank、minor、GB、jet又はCASは計算していない。

対象はSS_H_RESONANCEとNIL_H_RESONANCEの全P1・全九P2である。
SSはprojective rank二以下、nilはrank二以下又はactual quinticに一次因子を持つ。
従って両者のrank四・被約既約quintic部分は全て空。
旧G次元20/21を全点帰属へ読み替えた結論ではない。

## 1. 六shearは全域の実逆、全九P2を保持する

既独立receiptは元350×30/9、完全K1、実30×6 shearとそのfree-coordinate行を確認済みである。
SSの六shear座標は `(3t_a,3t_b,9t_c,-3s_a,-3s_b,-9s_c)`、
nilは `((3/2)s_a,(3/2)s_b,(9/2)s_c,(3/2)t_a,(3/2)t_b,(9/2)t_c)`。
係数は非零定数なので全parameterに実linear inverseがあり、zero chartを捨てていない。

有限変換 `a_i'=a_i+s_i*z+t_i*w` とその逆はdet一。
normal weightでP0は0/1/2、P1は1/2に移り、P2の全九成分は不変である。
weight一Euler regaugeの後も、P2は元P2にshear/P1の多項式を加えるtriangular変換で、
九変数についてのlinear partはidentity。
よってP2を零と仮定せず、全実解を残sourceへ移せる。

Eulerの変更は元projective Jacobiを保つ。
homogeneous quadratic bivectorはfull Eulerと可換なので、P0のEuler変更は
`E wedge [P0,P2]` というJ2を変えない。
従って旧J2のrank九は再gauge後も、**固定した実P1に解P2が存在すれば一意**という意味で使える。

## 2. SS全源のEuler展開

H=(a,b,-2c)、U=z∂z-w∂wとする。
`P0=H wedge U-5zw∂z wedge ∂w` から `E wedge U` を引く。
normal-normal部分は `2zw-5zw=-3zw`、c-normal部分は
`-3cz∂c wedge ∂z+3cw∂c wedge ∂w`。
従って正確に `-3Jac(czw)` となる。
残二列は `alpha0*z^2∂c wedge ∂w+alpha7*w^2∂c wedge ∂z` であり、

`Jac(-3czw-alpha0*z^3/3+alpha7*w^3/3)`

に一致する。これは全二parameterで厳密PoissonなのでP2=0が実解。
J2 injectiveから任意の元実解もP2=0に等しく、三変数Jacobianなのでprojective rank二以下。
実shearとEulerを戻してもprojective rankは変わらない。

## 3. Nil全四残parameterと正確なEuler再gauge

残parameter p,q,r,sについて、producerに表示された全P1から

`E wedge [z∂w+p*w∂a+q*w∂b-(r/3)*w∂c]`

を引くと、全sourceは

```
P=Jac(F)+(p∂a+q∂b) wedge W+P2,
F=(3/2)*c*z^2+(r/3)*z*w^2+(s/3)*w^3,
W=w*(-4c∂c+2z∂z+w∂w),
P2=A(z,w)∂a wedge ∂b+B(z,w)∂a wedge ∂c+C(z,w)∂b wedge ∂c.
```

と正確に一致する。
p列では `-bw` のab成分が消え、ac係数は-5+1=-4、az係数は1+1=2となり、aw係数w²が生じる。
q列も同様。r列は `(2r/3)zw∂c wedge ∂z-(r/3)w²∂c wedge ∂w` となりJac(rzw²/3)。
P0はJac(3cz²/2)、s列はJac(sw³/3)である。
この展開でa,b依存は全て消えるが、A,B,Cは全binary quadraticsのままである。

## 4. 存在する点だけでの一意性とactual stabilizer

GL2(a,b)はJac(F)、Wと全九P2空間を実際に保つ。
ここでは任意のP1にP2解が存在するとは仮定しない。
実解が一つ存在する点を固定すると、そのP0+P1を保つ群の像も同じJ2方程式の解なので、
点ごとの一意性からP2自体がその群に固定される。
chosen left inverseのequivariance又はgeneric compatibilityは要らない。

p=q=0ではGL2全体を使う。中心lambda*Iはab成分をlambda²、ac/bc成分をlambdaでscaleするので、
不変vectorは零、従ってP2=0。
この場合も三変数Jacobianでrank二以下。

(p,q)非零では、p非零chartの実matrix
`[[1/p,0],[-q/p,1]]`、p=0,q非零chartの `[[0,1/q],[1,0]]` がvectorを∂aへ送る。
逆はそれぞれ `[[p,0],[q,1]]`、`[[0,1],[q,0]]`。
両chartと零補集合を明示的に保持する。

∂aを固定する実stabilizer `[[1,t],[0,d]]` のdiagonal部分は
abとbc成分をd倍し、ac成分を固定する。
一意性からA=C=0。unipotent部分もこの残ac空間を保つ。
よって全実解はb-rowが零となる。
actual Euler-Pfaffianではb以外の全収縮項がb-rowを使って消え、残りは符号を除いて
`b*Pf(P restricted to(a,c,z,w))`。
したがって零又はproper linear factorを持つ五次式であり、被約既約quinticではない。

この全source空性を受理する。nilのP2を一般点で解けると仮定せず、schemeのradicalも用いない。
他のdependent N源への既essential三成分の帰属は別の既受理入力であり、本稿の新結論はH二源だけである。
global、成分数、manuscript、旧producer bytes、Gitは変更しない。
