# constant DOUBLE/TRIPLE の全 long25 親排除：最高symbolの共通 image-line

2026-11-08。PRIMARY E、担当 `/root/plane_collision_review`。
**新純証明producer、独立査読待ち。** 新point/operator/minor/GB/profile/局所変形jetの計算無し。
TRIPLEの接次元28は25次元親を排除しない。本稿は別の proper closed 全源条件を使う。

## 1. invariant marked-H と affine quadratic symbol

H⊂P4をPoisson invariant hyperplaneとする。
座標hをその方程式に取り、元homogeneous quadratic liftをΠとする。
H不変性はΠ#dh|HがEulerに比例することであり、h-rowはH上で
Πhi=λ*xi と書ける（λはlinear section）。
Euler∧linear場を引いてまずこのλを消すと、h-rowの全entryはhを割る。
さらにEuler∧linear場を引き、**全h-rowを零**に出来る。
これらはprojective tensorを変えず、h=1のintrinsic affine bracketを変えない。
同じ結果は affine chartから直接得られ、そのbracketはquadratic以下である。

h=1のaffine空間のtranslation vector space V_H はdimension四。
affine originを変えても最高degree二のsymbol

```
Q_H ∈ Sym²(V_H*) tensor Λ²V_H
```

は変わらない。affine frame変更でGL(V_H)共変である。
hの非零尺度とtensor尺度も以下のconditionを保つ。

共通image-line条件を

```
∃ l⊂V_H, dim l=1, Q_H ∈ Sym²(V_H*) tensor (l wedge V_H).     (L)
```

と定義する。すなわち全monomial coefficient bivectorsが同じconstant line lを含む。
zero coefficient、Q_H=0も許す。
固定Hでは、P(V_H)上のこのlinear incidenceを忘れる射がproperなので(L)はclosed。
H/frameを局所自明化すると上述gaugeと最高symbol抽出はlinear/regularであり、
origin変換による独立性から、(L)は全marked invariant-H incidence上でもclosedである。

## 2. 全long源は(L)を満たす

既受理必要十分sourceは、座標(X,Y,z,w)、r=Y−X²/2で

```
PXY=−F−(1+d)r, PXz=u, PXw=v,
PYz=Xu+r*g_w+H_w, PYw=Xv−r*g_z−H_z, Pzw=g,
d=div U, U=u∂z+v∂w,
deg F≤2, deg U≤1, deg g≤1, deg H≤3.
```

全Jacobi条件U(g)=dg,(U+1)H=Fgを保持する。
重要なのは、PXz,PXw,Pzwにdegree二項が**一切無い**こと。
従って全quadric coefficient entryはY-rowにしか無く、全域で

```
Q_H=∂Y wedge W2
```

となる。全F六係数、全H三次、gの定数/非定数/零、Uの全Jordan/translation/共鳴、
全十七源のzero/unit境界を捨てない。
rは原式を印刷する記号であり、rへの非線形変更をactual affine/PGL同値として使っていない。

従って全17long源のPGL像は、markされたmultiple hyperplane Hに対して(L)を満たす。
その全閉包でも次節により保持される。
各source閉包が既成分である行（μ1,μ1/3,非定数g,NEW25,M1/M11,B1–6等）では、
これはその**全閉成分**の必要条件となる。
short全源は既分母なしactual curveにより同じlong源の閉包に入るので同じ必要条件を持つ。

## 3. multiple-Hの proper closed flag と moving frames

rank四 locusでactualPfaffian quinticをDとする。

```
∃ H, H² divides D, H invariant, and condition(L)
```

を考える。H²|Dはprojective multiplication incidenceでclosed。
H不変性はtangencyのclosed linear condition。
(L)のline parameterもproperで、hyperplane parameterのdualP4もproperである。
従ってこの全存在条件はrank四parameter space上でclosed、PGL不変。
ゼロleading symbolも含むので、singular basisやperiod polesの境界を落とさない。

より具体的には、任意のlong源曲線でmarked Htはproperに延長し、
Ht→H0の後、DVR上のregular GL5 frameでHtを固定H0に出来る。
残るaffine frameが発散していても(L)自体はproper Grass/line incidenceだから保持される。
一つの固定coframeでのLaurent収束だけに仮定を狭めない。

## 4. TRIPLE の全monomial coefficientsには共通lineが無い

受理済みcanonicalのh=1 quadratic highest symbolは

```
Q=−xy ∂t∧∂y+(xr−3y²)∂t∧∂r+3x²∂y∧∂r.
```

xy,xr,x²の三monomial係数はそれぞれ

```
−∂t∧∂y, ∂t∧∂r, 3∂y∧∂r.
```

非零 decomposable bivector B に対しB∈l∧Vなら、lはそのimage二planeに入る。
上の三image planeの交差は

```
span(t,y) intersection span(t,r) intersection span(y,r)=0.
```

従って(L)は不可能。
TRIPLE D=h²*C3でC3はh²をleading coefficientとするr一次のprimitiveな絶対既約cubic。
したがってdouble hyperplaneは唯一H={h=0}であり、別のHで(L)を満たす逃げ道はない。
§3のclosed条件を満たさないから、**TRIPLEは全long-source閉包／対応する全long25成分の外**。
特に定数g・μ1/μ1/3の二25親候補を排除する。

## 5. DOUBLEにも同じ必要条件は失敗する

DOUBLEのquadratic symbolには

```
xy: −2∂t∧∂x+2∂y∧∂r
```

というrank四のconstant coefficient bivectorがある。
l∧Vの全bivectorsはrank二以下なので、この一つだけで(L)は不可能。
D=h²*C3も唯一のdouble-Hを持つので、DOUBLEも同じ全long閉包の外。
この補足は既D101による25以上親排除と整合するが、そのambient算術を再使用して証明していない。

## 6. 範囲と残余

本稿は、全source像閉包が成分そのものだと既証明されたlong25の全親排除である。
単にあるlong点を含む任意のより大きい別成分まで排除したとは言わない。
CA/CBやMPS成分が特定long境界点を含むことだけから、その全成分へ(L)を拡張しない。
CA/CB等の別factor閉条件、MPS T01/T39/T40、T36/T37/T38、Ga9などの別親監査は各既定理の範囲を保つ。
新成分性、全部分族のsmooth性、全E完了を主張しない。
