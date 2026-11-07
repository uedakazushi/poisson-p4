# normal CI(2,2)曲面上の全rank二restriction：既約quinticは存在しない

2026-11-07、主義務G、producer plane_collision_review。
**命題（producer証明）。** normal integral nondegenerate CI(2,2) S⊂P4上で
Poisson restrictionがrank二となるrank四Poisson構造は、実退化因子が既約被約quinticとなることはない。
restrictionの零divisor Rのreduced/irreducible/vertical条件は一切課さない。
non-normal CI22とprojected quartic scrollのrank二枝は本命題の外である。

## 1. 全sourceとsurface正規形の網羅性

`NORMAL_CI22_RANK_TWO_TANGENT_REDUCTION.md` §§1–3に全sourceの幾何を記す。
normal Gorenstein adjunctionとcodim二Hartogsにより、restrictionはOS(1)全五section。
constant four-vector UのJacobian bivector JUで全てlift出来る。
残りはKoszul/Euler完全列から全q1A+q2Bで、JUは両qのstrict Casimir、全UでPoissonかつJU²0。
全conormal条件は既受理のrankzeroの場合と同じN(M)二份である。

dimN≤1の全型ではΠ=JU+alpha C、alpha∈span(q1,q2)。
実D=E∧(2alpha JU∧C+alpha² C²)はalphaで割れるか零で、既約quinticに不可能。
これはPoisson条件を解く前の全sourceに対する因子である。

normal性とsize五のself-adjoint Jordan分割を使うと、dimN≥2の残りは次の全三surface orbit：

```
sem221:       q1=x0²+x1²+x2²+x3²+x4², q2=x2²+x3²+2x4²;
mixed21_11:  q1=2x0x1+x2²+x3²+x4²,   q2=x1²+x3²+x4²;
nil41:       q1=2x0x3+2x1x2+x4²,     q2=2x1x3+x2².
```

どの固有値もeigenspace dimension三以上ならrank二quadricを含み、非退化integral性に反する。
二blockの長さが双方二以上ならq1はそのkernel二空間上零、singular lineを生じnormalでない。
従ってc二の(22)+(1)、単一(3,2)はnon-normal補集合である。
残りc二はsem221とmixed21_11、単一固有値でdimN二となる型はnil41だけ。
pencil PGL2で三/二/一固有値を固定できるため、連続modulusをこの三表から落としていない。

各々の全N二generatorは順に

```
sem221:       R=∂0∧∂1,               T=∂2∧∂3;
mixed21_11:  R=∂0∧∂2,               T=∂3∧∂4;
nil41:       R=∂0∧∂4,               T=(3∂0∧∂3+∂1∧∂2)/4.
```

全sourceを `Π=Σui JUi+(a q1+b q2)R+(c q1+d q2)T` とする。
JU_iはmissing index iのconstant four-vectorを使う。全九係数を保持する。
二×二係数matrixのdet=ad−bcが零なら、二つのquadric係数は従属なので再びΠ=JU+alpha C。
したがってdet0の全補集合はquadratic factorで排除される。
以下det非零の全openだけを扱う。

## 2. sem221では全imageが固定四空間に入る

full original350Jacobiのexact rowspaceには次の式がある：

`ac=bc=0, a u0=a u1=0, d u2=d u3=0, d(a+b)=0`。

det非零ならc0、a,d非零。従ってu0=u1=u2=u3=0、b=−a。
残るJU4、R、Tのimageは全てW=span(e0,e1,e2,e3)に入る。
実D=E∧Π²はx4で割れる。従って全rank二openにも既約quinticはない。

## 3. mixed21_11でも全imageが固定四空間に入る

元Jacobiのexact rowspaceには

`ac=bc=0, a u0=a u2=0, d u3=d u4=0, d(a+b)−2b u2=0`

がある。det非零から同じくc0、a,d非零、u0=u2=u3=u4=0、b=−a。
残るJU1、R、Tは全てW=span(e0,e2,e3,e4)に入るので、実Dはx1で割れる。

## 4. nil41はdet非零openが空

元Jacobi rowspaceに `c²=0, 2ad−5bc=0` がある。
複素数上の被約支持でc0、ad0となり、det=ad−bc0。
従ってdet非零openに点はない。cのscheme nilpotenceを零と置換した全scheme同型は主張しない。
この主張はPoisson点の被約分類に十分である。

## 5. regular member無しのnormal補集合も空

all-symmetric-pencil canonical theoremの全coverは既独立確認済みである。
M1+regular2のintegral modelはsingular lineを持つ。
common radical一ではbaseがelliptic quartic C⊂P3であり、C singularなら再びsingular line。
従ってnormalとして残るのはsmooth Cのconeだけ。
全modulusをq1=Σxi²,q2=Σλixi²（四λ相異）で表す。

全conormal sourceはΠ=JU+∂v∧(q1V+q2W)、V,W constant base vectors。
Uのbase-four-vector係数をcとする。c0なら全Πが∂vを共有しrank二以下。
c非零なら、両qのCasimir性とJU全compatible性により、元Jacobiは

`Ebase∧(q1 Lie_V Jbase+q2 Lie_W Jbase)=0`

を必要とする。missing dのcomponentのxi⁴係数は、非零共通scalarを除き
`(λd−λi)(Vd+λi Wd)`、i≠d。
異なる二つのiを取るとVd=Wd0、全dでV=W0となる。
するとΠ=JUで全rank二以下、実D0。全elliptic modulus/vertex係数を保持している。

## 6. 証明書とscope

全三型のraw25/Euler24はcombined rank49、全conormal spaceは二generator。
全九sourceの元150係数、350Jacobi、実126Pf、equivalent quadratic rowspaceを
`CLASSIFICATION_NORMAL_CI22_RANK_TWO_SOURCE_20261107.json` に保存した（0.147秒）。
補完checkerは§§2–4の各必要式を**元350行のexact functional**として復元し、
sem/mixedの全縮約三parameter Poisson源のactualDにx4/x1因子を直接確認する。
これらは例点でもjetでもなく、全元sourceの恒等式である。

normal性に依存するOS(1)延長をnon-normal曲面へ広げない。
E帰属・新成分数・全Poisson schemeの非被約構造は本結果から数えない。
一次文献のcanonical-pencil入力以外の新幾何は上記で証明したが、新規性は未監査。
