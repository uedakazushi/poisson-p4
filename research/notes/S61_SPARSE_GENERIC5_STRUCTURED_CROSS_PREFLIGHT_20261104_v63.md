# cross126の構造的逆に向けた実第一列の縮約と残る84列

2026-11-04。新しいfield Gaussian、Gram係数展開、数値β、素数、GBを使わず、実adapterの一列を簡約し、小さいpolar blockのsupportを測る。全left inverseの係数は未計算である。

## 1. 実E第一列は75項の既知差分になる

K=QQ(H,L,s0,q0,t0)、d=Q52-centerとする。保存された45座標を直接比較するとQ51=center-d/2である。現行q pivot行は[3,6,13]で、厳密に

```
d_pivot=(H*q0*t0/144,0,0).
```

受理済みW=P*[NT0,NT1,d]、E=W*adj(Wpivot)/Gammaに対し、普遍逆行列恒等式からE*d_pivot=d。従って

```
E[:,0]=144*d/(H*q0*t0).
```

d自体は45成分中41成分が非零、合計75項のQQ五変数多項式である。Eのこの一列に、以前の大きいCramer/NT/Gammaを再展開する必要はない。反転するのはparameter単項式だけである。

## 2. 134×42の小さいraw polarと、quotientの実cost

raw F134のq二次係数をQ_ijとすると、d方向と現在のstd_nonpivot法線n_kとのpolarは、Q_ij*d_jを足すだけで得られる。対角は二倍する。全入力係数を計算し直さず、monomial supportの和集合を取るpreflightで次を得る。

|block|全support上界|一成分上界|parameter次数上界|
|---|---:|---:|---:|
|d方向raw polar|13506|14|6|
|Q51でのnormal Jacobian|14108|14|6|
|Q52でのnormal Jacobian|14081|10|6|

13506は相殺前の上界であり、新しい全QQ係数を構成した数ではない。共通numeric分母を使う安全な係数上界は分子89bit・分母44bitである。

しかしraw polarを残差crossと混同してはいけない。実C_d,resはraw polarから既存五本のselected polarをlambdaで引いたもので、lambdaには実五汎関数selectorが必要である。このactual K5 selector入力は**692953項**で、150000項のfield-stage guardを超える。外部回路として参照することはできても、13506だけを全K5入力costと呼ばない。

同じformal QQ[T]曲線に潰した場合だけはselector55594項＋raw polar上界13506=69100項、raw polar次数20以下、47minor次数1395以下になる。この数は五変数係数costと別物である。rank47自体は既受理のcross126から従うので、既知rankの再確認だけを目的とする新formal rank jobは提案しない。固定小minorや実係数逆に結び付く新しいgateが先に必要である。

## 3. 正確な42＋84分割recipeと、未計算gate

既存actual X130×126を[X0|X1|X2]と並べると、各block幅は42である。前節からX0=C_d,res/(H*q0*t0/144)。この42列の全source rootを別証拠に保存する。rankX=126からX0は単射だが、具体的な42minorやscalar逆係数は計算していない。

形式的にはZ=X0、B=[X1|X2]と置き、parameter-onlyな

```
LZ=(Zt Z)^(-1)Zt,
Y=(I-Z LZ)B, LY=(Yt Y)^(-1)Yt
```

によりyの逆像を`(LZ*(y-B*LY*y),LY*y)`と書ける。Kはformally real、Xは単射なので二つのGramは可逆である。ただしこれは126を42と84へ分ける有限演算recipeであり、**構造的scalar逆を発見したという結論ではない**。射影が残る84列を混合し、二つのGramの係数増大上界も未証明だからである。従来の126 Gramの代わりにこのrecipeを実行したとも、総計算量が必ず下がるとも主張しない。

## 4. spectral分解・五点補間・homotopyだけでは足りない理由

一つのtoric Poissonテンソルを対角化すれば、そのPoisson differentialにはmonomial weightごとのKoszul記述がある。しかしcrossは三つの独立なplane方向と42法線のtensor積を写す演算子で、単一テンソルのdifferentialではない。三つのtoric点が同じlogcanonical対角frameを共有すれば任意の線形結合もPoissonになり、既知planeとの交わりが孤立した被約三点であることに反する。fixed betaのS61再構成はqについてaffineなので、この反証をq平面にも適用できる。

五点代数の冪等元はplane上の関数を分解するが、normal42や方程式空間に対する作用は自動的には与えない。R_h(0)をその作用と置くと、それら自体が未知のcross左逆に依存し循環する。Clifford C/Dの一方の固有frameを取っても、他方との可換性やcross上の共通weight分解はまだ証明されていない。

有限分離拡大L/K上で**全恒等式付き**左逆Yを得た場合はTr(Y)/[L:K]がK上の左逆になる。frameの行列変更も元のXへ戻してからtraceする必要がある。traceは降下のgateを解決するが、L上の小block分解と係数formulaを代わりに与えるものではない。

また定数symplectic部分Jのde Rham contractionは係数の空間次数を一つ上げる。高いJacobi出力へ適用して得る二重ベクトルcochainは、固定S61 q45の許された二次・三次再構成空間に一般には戻らない。既受理inverse-form Bianchi/Koszul恒等式も、このcoefficient inverseを与えない。unknown空間変数のPfaffianをparameter単元として反転してはいけない。

## 5. 次の有限目標

有用な新入力はE第一列の75項formulaとraw d polarの13506項costである。次の目標は、五点標準模型への実adapterまたは低係数の五次元plane-image frameを束縛し、第一42列のleft inverseをactual parameter係数で作ること。その前に692953項selectorを再度field Gaussianへ入れない。二つのtoric方向との84列結合、normal Gram126、残余の空性は未解決である。

source/support確認は0.3713秒で完了した。新しい134×42係数多項式、42minor、field逆行列は構成していない。既存Xの第一42列全5460個のsource rootを保存し、第一列縮約の75項入力と全45座標差分を束縛した。
