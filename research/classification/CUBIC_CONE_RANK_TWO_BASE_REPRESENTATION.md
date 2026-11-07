# twisted cubic cone rank二：全base支持を五whole moduleへ

2026-11-07、PRIMARY G、root producer。全rank二Gの結論は最後のmodule確認後だけに採録する。
actualD既約被約quintic、twisted cubic cone S=(x0x2−x1²,x0x3−x1x2,x1x3−x2²)、vertex c=x4。
S上一般rank二ならnormal multiplicity二・abelian conormalが必要。全vertex/全Rを保持する。

## 全二十七源と元全三百五十式

全150 tangency及びintrinsic abelian firstnormal equationsは975×150 rank99、完全kernel51、
Euler24を除きcanonical divergence-zero源27。初回dimension guard27では止まり、G成功としなかった。
この完全源はliteralに

    Π=P(x0,x1,x2,x3)+(c div(P)+V(x)) wedge ∂c

となる。Pは全九次元、Vは全十八次元のdivergence-free quadratic field。
base indices [0,1,2,3,4,5,6,7,16]、他全18がvertex source。
全V-V Jacobiは恒等零、元350のpurebase二次及びbase×Vは異なるspatial coefficientに現れる。
したがって全方程式は全九変数の十二quadricsと、全V十八の350行homogeneous linear moduleへ分離する。
元全126actualD係数もbase二次＋base×Vとして保存した。
`CLASSIFICATION_CUBIC_CONE_RANK_TWO_BASE_VERTEX_20261107.json` は0.536秒で全恒等式を確認した。

## 実SL2と三つのbinary quadratics

coneのactual SL2は四base座標のSym3表示、cを固定する。全150係数でのLie h/e/f作用を計算し、
九sourceの正しいthree Sym2表示を同定した。t0,…,t8は上の順のbase coordinates。

    U=(3t0−t5, 3t4−2t6, 3t2+3t3−t7),
    V=(2t0−t5−t8, −4(t1−2t4+t6)/3, (3t2+4t3−2t7)/3),
    W=(t0,t1,t2+t3).

各三つ組(a,b,c)はbinary quadratic a s²+b st+c t²の係数である。
全9×9変換と逆行列、各blockのh=diag(2,0,−2)、f:(a,b,c)↦(0,2a,b)、e:(a,b,c)↦(b,2c,0)を確認した。
十二Jacobi rowspaceはUVの五coefficientsを含み、全元十二行からの明示functionalを保存した。
C[s,t]がdomainなので、全reduced supportでU=0又はV=0である。
ここでtarget coordinateやunknown source coefficientを反転していない。

残りの全十二式を各whole branchへ代入し、両方向rowspace同値を得る。

* V=0：UとWの三2×2 minorsが零。
* U=0：VとWの三2×2 minors及びdisc(V)=v1²−4v0v2が零。

したがって全支持は、V0でrank[U,W]≤1、又はU0でrank[V,W]≤1かつVがnil quadraticである。
`CLASSIFICATION_CUBIC_CONE_RANK_TWO_BASE_SUPPORT_20261107.json` は0.255秒で全変換/全式を確認した。

## 五canonical familiesと全補集合

U≠0の枝ではW=λU、λは全complex line。SL2とtensor scalarによりUは相異二根q_ss=st又はq_nil=s²。
これはactual ambient transformationであり、vertex十八源も同時に可逆に変わる。
U=0,V≠0ならVは平方でW=λV、同じactual操作でV=s²へ。
U=V=0ならWは相異二根、平方又は零。W零ではbaseP零、Π=Vfield wedge ∂cのrankは二以下。
これで全rank四sourceは三whole lines及び二fixed canonical baseへ尽くされる。

| case | canonical base | 全base stratum dimension | fixed λでのvertex kernel kに対する全PGL像上界 |
|---|---|---:|---:|
| A_ss | V0,U=st,W=λst |4|19+k|
| A_nil | V0,U=s²,W=λs² |3|18+k|
| B_nil | U0,V=s²,W=λs² |3|18+k|
| W_ss | U=V0,W=st |3|18+k|
| W_nil | U=V0,W=s² |2|17+k|

coneのPGL stabilizerは8、PGL orbitは16、tensor projectivizationで一を引くので上表となる。
λ有限wallではbase dimensionが一減るため上界も一減る。全zero complementsは上の別caseへ保持した。
二項行列のrankone記述をdense openだけへ局所化していない。

次の一finitejobは各全350×18 moduleのgeneric rank/完全核/非零minorを求め、minorとkernel分母の全irreducible factorsを処理する。
各wallは元matrixを戻す。dimensionが十分ならそれだけでG≤24、そうでなければ全actualDにproper共通因子があるかを確認する。
未検査wall又は非閉kernelは未解決と明記する。新GB/jet/追加point列を用いない。
