# Q1/Q2の二次functionalはrank44：T1=T2=0

2026-11-04。既存の二つの57-systemの全二次データだけを再利用し、
残る二次classを判定した。独立監査用の全係数・minorを保存した。
元B18原点のfull D depth/nilやT0はまだ決定していない。

## 1. 判定する最小関係空間

受理済みT3消滅の構造定理は、安定kernel Kinfの最小generator数が44、
その一次部分のrankが43で、残るclassがN moduleの二次関係に対応することを与える。
Q1/Q2は同じhonest9 frameへ正確に持ち上げ済みで、両方の一次spanも43。

同R0上のNのCohen座標を(z1,...,z8,v)、N generatorをn0=1と七つのodd generatorとする。
関係r∈Kinfをactual B_NでN generatorに送ったとき、そのn0係数について

```
phi(r) = [v²] (n0係数をz1=...=z8=0へ制限).
```

Nの全relationは、35本の線形関係と例外的な
`(v²-nu*z1²)n0`の一つで生成される。
全relationのn0係数はm_R0²に入るので、phiはm_R0*Kinfを消す。
43本のcanonicalな線形liftとnil curveの8-CI由来関係にphiは0、
例外的な二次liftにphiは1となる。
従って一次43座標にphiを加えると、Kinf/mKinfの44次元全体を検出する。
このdescentと最小性はT3構造定理および同R0 N module36関係を使い、
一般のベクトルのorderだけからmembershipを推測しない。

## 2. 必要なのはphysical x61線上の係数だけ

actual R8はphysical rows
108,109,112,113,117,118,119,27であり、rho61=a61(R8)+v。
新z8はR8のunit座標、a61(0)=0なので、z8=0は正確に

```
physical R0の八Cohen座標=0,  x61=v
```

という一本の線である。一次近似だけの同一視ではない。

B_Nのfirst rowは、honest generator jに対してa_h(R8)+c_1j(R8)*v。
従ってこの線の一次係数はc_1j(0)だけであり、保存済みの値は

```
(0,0,0,0,0,0,0,1,-1/2).
```

これはCA evenのNへのbase changeも含めた、全honest generator像の係数である。
raw physical odd係数に置き換えていない。
old relationの定数は全て0、B_Nのfirst row定数は(1,0,...,0)なので

```
phi(r) = old_E0_second_on_line
          + sum_(j=1..8) c_1j(0)*old_r_j_first_on_line.
```

nil補正行kappaは不要である。

## 3. 二次の座標変換を正確に含める

CA/CBそれぞれの保存済みsource二次jetは、各自のF9座標Yで記述されている。
共通physical x(Y)=L*Y+Q(Y)+O(Y³)から、x61線の逆像を

```
Y=y1*v+y2*v²+O(v³),
L*y1=e61,      L*y2=-Q(y1)
```

として正確に求めた。
全source係数f(Y)は、保存済みの一次・二次係数をこの式へ代入してv二次まで戻す。
二経路のphysical x9が、定数0、一次e61、二次0に厳密一致することを確認した。
二次逆変換を捨てた線形代入ではない。
source IFTの再実行、三次source jet、whole honest二次matrixは作っていない。

## 4. 未計算の正則商の微分を消したことを検証する

旧coordinate商784個の一次係数は、FULL9で既に正確に求めたものを使う。
それらの二次係数と四つのg商の二次係数には、独立なweight2形式変数thetaを割り当てた。
四つのg商の一次係数にも独立なweight1形式変数gammaを保持した。
これらを0と仮定していない。

raw E0を全Q1/Q2列で二次まで計算し、unit Schurの式

```
new_E0 = raw_E0 - pivot_even * old_component5 / unit
```

を適用した。pivot_evenの原点値は0、unit原点値は10/3。
一般の二次積公式を使い、old_component5とunitには一次だけを使った。
実際にはpivot_evenの線上の一次部分も正確に0だった。

全phiにおいて、全thetaと全gammaを含むmonomialの係数が**全て厳密に相殺した**。
従って全1548/1961個のphiは正確な有理数となり、未知微分を捨てた上界ではない。
全形式変数対応と各列の相殺後の式を保存した。

## 5. 44次minorとNakayama

全一次90係数にphiを付加した91成分column行列を作った。
Q1、Q2ともrankは **44**。
選択した44×44有理minorのdeterminantは、両方とも

```
-390625/279936 != 0.
```

したがってQ1,Q2の像はKinf/mKinfの全44次元を張る。
有限module Kinf/QnへNakayamaを適用して

```
Kinf=Q1,      Kinf=Q2,
T1=0,        T2=0.
```

が従う。T3=0は既存の独立結果であり、今回の計算で残る正次数二層も消える。
ただしT0は別moduleである。w-torsionがT0だけになっても、その有限台torsionを
まだ判定していないため、full Dのdepth19/20やnilmoduleの全分類は主張しない。

## 保存・計算量

Source `code/ambient_B18_T1_T2_quadratic_phi_20261104_v63.py`。
`reports/AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63/quadratic-phi-certificate.json`に
各族の正確な二次座標逆変換とsource係数線jet、全形式微分の対応、全列のphi、
全augmented matrix、選択44次minorを保存した。
一回3.168435秒、最大RSS459400KiB。60秒/2GB/32MB内。
新しいsource IFT/jet、whole honest二次matrix、GB、点、p逆元は無し。
