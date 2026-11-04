# T1・T2二次functionalの独立監査：両方零

2026-11-04。作者producerを再実行せず、既に独立検証した二つの57-sourceの保存二次係数を使う。
physical x61線への二次逆座標変換、全未知微分の相殺、44次minorを別の演算表現で再構築した。
actual T1=T2=0を受理する。元B18のT0・全D depth19/20・全nil構造は未決定。

## 1. 全次数の判定空間

独立受理のactual N同R0結合とliteral nilcurve extension、および
`AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63`の44生成証明を使う。
安定kernel JのJ/mJは44次元で、一次像が43次元、一次kernelは例外関係r0のclassのみ。

N-unit generator n0の係数は、Jの全generatorでm²に入る。
従ってその二次係数functionalはmJを消し、J/mJへ降りる。

```
phi(r)=[v²] (N-unit係数をz1=...=z8=0へ制限).
```

canonicalな43一次generatorにはphi=0、r0にはphi=1。
一次43座標とphiの組はJ/mJを同型に検出する。
43generatorの別の選択はphiへの一次座標のshearを起こし得るが、combined rankを変えない。
またB_Nのlift選択の差はN関係であり、そのunit係数はm²なので、B_Nのunit-row firstjetも不変。

## 2. 同じhonest9 frameとexact Q1/Q2入力

全raw Q1/Q2 trace pairを凍結済みQ3 pairへ照合し、各pairのliftが同じhonest9列に
一致することを再確認した。元36 Koszulを全て保持する。
unit10/3のpivot pairはQ1 index1243、Q2 index1353にも入り、同じ一行消去が正当である。
全列indexとmultiplicityを、unit消去後の元Q3 indexへ直接戻して照合した。
関係数は1548、1961であり、任意の関係を欠いたspecializationではない。

## 3. Physical x61線への二次変換

actual R8はphysical rows108,109,112,113,117,118,119,27そのものであり、
新z8はそのunit座標である。z8=0ならR8=0、a61(0)=0なので、old physical x61=vが正確に成り立つ。

各族の保存済みsource(Y)二次係数から、x(Y)=L Y+Q(Y)として

```
L*y1=e61,   L*y2=-Q(y1),   Y=y1*v+y2*v²+O(v³)
```

を独立Gaussで戻した。全57保存係数をこの式でcontractし、共通physical x9が
定数0、一次e61、二次0になることを両族で検証した。
作者保存の全contracted coefficientsとも一致した。
新しいIFT/Newton解法やsource二次以上の計算は行っていない。

## 4. FULL-HのB_N補正を独立に導く

96/160系producerを再実行したり、保存したc1値をそのまま仮定したりしない。
actual old honest generatorのrow jは

```
H_j=physical_rho_j-a_CA,j(R0).
```

Nの八Cohen parameterはeven、残るrho61だけがa61+b61*tauである。
従ってそのN odd coefficientの一次勾配は正確に

```
d(b_Hj)=d(b_N,j)-(partial_rho61 a_CA,j)(0)*d(b61).
```

高次のa_CA Taylor項は、このodd coefficientではorder2以上へしか寄与しない。
CA even derivativeは、上で再構築したexact source lineから得た。
Nの全odd勾配は独立受理済みN bindingから再利用し、そのunit8 matrixを独立に逆行列化した。
このfull-H勾配をz座標へ移すと

```
(c1j(0))=(0,0,0,0,0,0,0,1,-1/2)
```

となり、作者値に一致する。CA even subtractionを落としたraw physical odd係数ではない。

## 5. 別表現による未知係数の全相殺

作者はt/gamma/thetaを一つのweighted monomial環で扱う。
独立側は、ordinary t-jetの各係数をgamma/thetaの疎有理多項式とする別の表現を使った。
旧coordinate商784個の二次係数とg商4個の二次係数に788個の独立thetaを保持し、
g商4個の一次係数に四つの独立gammaを保持した。
これらの値を0へ置いていない。

raw evenとpivot evenを全列で二次まで評価した。Schurは一般の積公式

```
new_E0_2=raw_E0_2-p2*c0/u0-p1*c1/u0+p1*u1*c0/u0²
```

を使い、p0=0、u0=10/3、さらにp1=0をexactに確認した。
component5とunitは一次しか使わない。
それらの計算に新Hadamard792個の二次係数を置いたように見えても、二次まで必要な
raw even/pivot evenは旧full186 DAG内にあり、新792 nodeを一つも参照しないことを検証した。
新792の二次を任意値にした推論ではない。

全1548/1961列のphiにおいて、全theta、全gamma、gammaの全積の係数が厳密に相殺した。
従って各phiは正確な有理数となった。作者の全1961個のdistinct列の値と一致した。

## 6. Minorと全次数の結論

一次90係数は、私のT3構造監査で独立再構築した全2069列から再利用した。
各phiを付加した全91係数columnを作者matrixと照合し、独立疎Gaussでrank44を確認した。
保存された44次submatrixをentryごとに再構築し、独立有理Gaussで両方

```
det=-390625/279936 != 0
```

を確認した。J/mJのdimensionは44なので、両Q1/Q2が全てを張る。
NakayamaによりQ1=Q2=J、したがって**actual T1=T2=0**。
T3=0とdegree4安定化を合わせると、残るw-torsionはdegree0のみでwにannihilateされる。
そのT0の零・非零や有限socleをまだ判定していないため、全Dのdepth19/20は未決定である。

Code `code/ambient_B18_T1_T2_quadratic_phi_review_20261104_v63.py`。
Report/folder `AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW`。
二次coordinate逆変換、full-H correctionの独立導出、全phi、二つの44次minorを保存する。
