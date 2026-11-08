# 残m3二源の全Ga二対称性と一般orbit二十二・商一

2026-11-08、PRIMARY E、`/root/quartic_sources`。純手恒等式。新しいmatrix、rank、point、operator又は算術実行は無い。旧一点M47だけから一般orbitを推測せず、全source上の上界を先に証明する。

## 1. 固定する全sourceとstrict lift

`TWO_LINES_PLANE_CUBIC_THREE_MAGNETIC_QUOTIENTS.md` §§2–3の正規化でbase(h,x,y)、fiber(z,w)、h-Casimir quadratic liftを使う。全六磁気R2を残し、

```
P=U wedge partial_z+V wedge partial_w+(-d*w+R2)*partial_z wedge partial_w.
```

nodeではU=(0,2x(x-h),y(3x-h))、V=(0,-2xy,2x(h-x)-3y²)、d=x-h。cuspではU=(0,2x(x+h),3xy)、V=(0,2xy,3(y²-(h+x)²))、d=xである。

旧証明書 `CLASSIFICATION_TWO_LINES_PLANE_CUBIC_20261107.json` の両caseは全R6について十strict Jacobi零を保存している。したがって、以下でdivergenceのBV恒等式を使うのは**このstrict h-Casimir lift**である。任意のEuler代表のdivergenceをprojective stabilizerと呼ぶ誤りではない。

## 2. 全R6に依らない二つのlinear Poisson field

規約 `Z^j=sum_i partial_i Pij` で全五変数divergenceを取る。strict JacobiとBV恒等式から `[Z,P]=0`。また

`X=P sharp dlog(x/h)`

はclosed logarithmic one-formのHamiltonian fieldなので、h*x非零開でPoisson vector fieldである。xはPoisson-normal、hはCasimirだからXはglobal linear fieldとなり、全polynomial恒等式として延長する。

nodeでは直接微分により

```
Z=(8x-4h)*partial_z-8y*partial_w,
X=-2(x-h)*partial_z+2y*partial_w,
Z+4X=4h*partial_z.                                  (1)
```

cuspでは

```
Z=(8x+2h)*partial_z+8y*partial_w,
X=-2(x+h)*partial_z-2y*partial_w,
Z+4X=-6h*partial_z.                                (2)
```

例えばnodeのZ_zは `(4x-2h)+(3x-h)+(x-h)=8x-4h`、Z_wは−2y−6y。cuspのZ_zは(4x+2h)+3x+x、Z_wは2y+6yである。R2はfiberに依らないため全係数がこの微分へ寄与しない。Xの式もR2と無関係で、x-rowをxで割った全式そのもの。

両fieldはbaseを固定し、fiberへbaseのlinear式を加える。したがってcommuteし、その各行列の二乗及び互いの積は零である。有限actual flowは `I+sZ+tX`、det一、逆はparameter符号反転となる。両者は実際のGa² subgroupを与える。

独立性も全域で明示できる。nodeでaZ+bX=0ならw-rowからb=4a、z-rowから4a*h=0となりa=b=0。cuspもb=4aの後−6a*h=0となる。Euler scalar moduloでも独立である。両fieldのh,x,y成分は零なので、その定数線形結合がscalar Eulerならh-rowからscalar零となり、上の独立性へ戻る。

## 3. 全source上界と既nonzero minorによる一般下界

PGL5の次元は24、全sourceのstabilizerは上のGa²を含む。従って**全六磁気sourceの全点でactual projective orbit次元は高々22**。actual GL/PGL変更とtensor scalarを戻してもこの上界は不変である。作用matrixのrank条件は閉なので、そのPGL飽和閉包上にも対応する上界が延びる。

旧fullambient証明書の `cases` 内、

```
node_flex_tangent_through_flex,
cusp_flex_secant_through_flex
```

では各々 `ranks_QQ={D:101,vertical:24,orbit_scale:47,full_source:48}`、orbit_scaleの47minorと全sourceの48minorが保存されている。このminorの非零開は元A6で稠密である。Euler24とprojective scale一を除いたactual orbit下界は47-24-1=22、全source飽和像下界は48-24-1=23となる。

全source上界22とこの開下界を合成して、**各Siの一般actual PGL5 orbitは正確に22次元**。全像次元23の上界は既 `TWO_LINES_PLANE_CUBIC_KNOWN_PARENT_REDUCTION.md` §1のessential base Gr(2,5)六＋marked curve orbit八＋fixedD source九による全parameter数え上げであり、同じ48minorの下界と一致する。

従って両Siはdimension23、一般reduced Rosenlicht quotientはdimension一。具体的なquotient関数又は全lambdaの同値関係は分類していない。orbit22が全境界で一定とも主張しない。

一点D101からgeneric Poisson-scheme transverse ringを推測しない。一般lambdaでDが102へ上がる可能性等は本稿では計算していない。被約componenthoodと一般orbit/商に必要なのは上の全対称性・既minor・全source dimensionだけであり、新しいrank点やnil指数は不要である。
