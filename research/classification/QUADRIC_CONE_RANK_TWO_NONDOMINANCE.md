# quadric cone上の全rank二枝：PGL像二十四上界

2026-11-07。主義務G、producer plane_collision_review。
**命題（producer証明、別作者独立review前）。** 実退化因子Dが既約被約quinticで、
その非正規曲面supportがquadric cone Sであり、πのS上の一般階数が二である全Poisson構造の
PGL5飽和像はdimension二十四以下である。
S上のanticanonical零divisor Rの既約性・被約性・normal-cone型は仮定しない。
全成分のE帰属、存在・不存在、schemeの非被約構造はこの命題の範囲外である。

## 1. 全sourceと補集合

座標(u0,u1,u2,v,h)、q=u0u2−u1²、S=(h=q=0)。
`QUADRIC_CONE_RANK_TWO_TANGENT_SOURCE.md` §§1–3の全tangent Euler liftにより

`Π=J(q,F)+hL+qA`

と表せる。Fは全H0(OS(2))九係数、Lは全五十係数、Aは全十係数。
raw69、Euler冗長五、projective64であり、全150係数のrankとEuler24のcombined rank88も確認した。
conormal bracket零の全十一線形条件を保持する。
A(h,ui)0、A(h,v)=aとなり、a0ならHがactualDを割るため排除する。
a非零の全actual inverse shear/scalingsによりA(h,v)=1、A(ui,v)=0とし、
残るCをsemisimple ∂0∧∂2、nilpotent ∂0∧∂1、zeroの全三型へ直す。
全F sourceは同時に可逆に移り、F0も安全な過大評価として保持する。

各Cで全五十九unknownのmodh original350Jacobiとconormal十一式を解き、
raw affine dimension38、Euler quotient後33の**完全**kernelを得た。
全kernel方向、inverse free-coordinates、Euler minor、全150tensors、元350Jacobi、実126Pf係数を
`CLASSIFICATION_QUADRIC_CONE_RANK_TWO_PENCIL_20261107.json` に保存した。

## 2. 誤った入口を除き、元式からOreを得る

初回はrankzeroと同じv非依存Bを入口で仮定し、全modh source上のassertに失敗した。
その失敗source/job/runtimeは不変保存し、成功とは扱わない。
正しい全33源は

`P=B0(u)+vD+(Y(u)+vT(u))∧∂v`

である。Dはconstant bivector。全五つのHamiltonian-pencil Jacobi blockを
元350projective係数へ全parameterで比較して同値を確認した。
特に`Jac(B0)+Y∧D=0`を`Jac(B0)=0`へ勝手に強めていない。

`CLASSIFICATION_QUADRIC_CONE_RANK_TWO_HIGHEST_SYMBOL_20261107.json` の最高spatial係数から、
F_v=ell0u0+ell1u1+ell2u2はCsemでell0=ell2=0、Cnilでell0=ell1=0と従う。
C0では全三係数が残るがDは恒等零である。
従って**元Jacobiの帰結として**全三型でD0が得られ、全Ore条件

`Jac(B)=0, Lie_T B=0, Lie_Y B=Y∧T`

へ正当に戻れる。CsemではF_v=κu1、Cnilではκu2とする。

## 3. 全graded fiberを三つのlinear operatorで上から押さえる

全sourceをB2+B1+B0、T1+T0、Y3+Y2+Y1+Y0と分ける。Y3=−u q固定。
baseはκとradial trace s。全他の自由係数は
weight1十四（B1,Y2,T0）、weight2十二（B0,Y1）、weight3三（Y0）である。
元式の各係数はparameter weightでhomogeneousであり、各stageの新unknownに対してlinear、
そのdiagonal matrixの係数はbaseだけに依存する。
従って全Poisson fiberのdimensionは三matrixのnullity和以下である。
低degreeのquadratic compatibilityはfiberを小さくするだけで、この上界を増やさない。
全元rowとの対応とweight identityを確認し、全matrixを保存した。

| model | 三operatorの大きさ | generic nullities | PGL像上界 |
| --- | --- | --- | ---: |
| Csem | 38×14,21×12,7×3 | 3,1,0 | 24 |
| Cnil | 35×14,22×12,8×3 | 3,1,0 | 23 |
| C0,F_v sem | 33×14,21×12,7×3 | 3,1,0 | 23 |
| C0,F_v nil | 33×14,22×12,8×3 | 3,1,0 | 22 |

C0,F_v0ではfull126Pfのh0係数は恒等零なのでactualHfactorとなる。
非零F_vは残るSO(q)×scalarの全二軌道u1/u2へ移る。
このactual normalizationとinverse、全三/二次元orbitの回復は
`QUADRIC_CONE_RANK_TWO_TRIANGULAR_BOUND.md` §3を参照。

## 4. 全共鳴とnormalizer次元を戻す

Csemのselected-minor productの異なる全因子は次の十二lineである。

```
2κ+1, 2κ−1, κ, κ+1, κ−s, κ+s,
κ+s+1, κ+s+3, 3κ−s, κ−s+1, κ+s+2, 3κ−s+1.
```

Cnilでは七line：

`2κ−3, 2κ+3, κ−1, s, 4κ−1, κ, 2κ−1`。

全lineをQQ(tau)で解き、そのgeneric minorの全複素根をexactに戻した。
Csemは計七十二root certificates、Cnilは十三（重複する交点も意図的に保持）。
C0,F_vsemの全根s=−3,−1,1,3と、F_vnilの全根s=0も保持した。
全matrixの完全kernelとminorを保存し、generic minorの非零領域だけで結論していない。
base dimensionはgeneric二/line一/point零、C0では一/零として加算する。

AutS×tensor scalarは十三次元。Aのsem/nil/zero strataは七/六/四次元。
Asem/nil固定後のfull source dimension dに対する像上界は18+d/17+dである。
C0でさらにF_vのsem/nil orbitを固定した後も、orbit dimensions三/二を戻して18+d/17+dとなる。
全framed sourceのfree normalizer fiberを引いており、既に消費したnormalizationを二度引いていない。

この全coverで唯一Cnil,κ=s=0だけがcoarse image25となる。他の全strataは二十四以下。
算術は `CLASSIFICATION_QUADRIC_CONE_RANK_TWO_TRIANGULAR_BOUND_20261107.json`、
成功を偽装せず`all_rank2_quadric_cone_G24_closed=false`のまま不変保存する。

## 5. 唯一の全共鳴源は元compatibilityで二十四へ落ちる

Cnil,κ=s=0の三nullitiesは(5,3,0)。全top十四係数を完全kernelのp0,…,p4で表す。
元weight二の全二十四式はmiddle十二係数にlinear、rank九。
そのleft-annihilatorを取り、保存済みpencil-row indices13,59,62,75に
それぞれ7/8,−1/8,1/8,1/8を掛けるとmiddle係数は全て消え、正確に

`p1*(p0+p2/2)=0`

を得る。このidentityと全元matrix・全quadratic右辺を
`CLASSIFICATION_QUADRIC_CONE_RANK_TWO_RESONANT_COMPLETION_20261107.json` に保存した。
これは全五自由係数上の非零元条件であり、top supportはdimension四以下。
middle fiberは三以下、bottom fiberは零なので全sourceは七以下、PGL像は17+7=24以下となる。
free coefficientを一つも反転せず、二つのlinear factorsの零枝も双方保持する。

以上で全A/F_v/base共鳴と全restriction零divisorRを尽くし、命題が従う。
新規性・既存文献との同定は未監査。これはG側の全像dimension上界であり、
未知点のE帰属や新成分の数はこの結果から決定しない。
