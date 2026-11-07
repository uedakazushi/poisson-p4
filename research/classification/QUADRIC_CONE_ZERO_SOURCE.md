# quadric cone rank零枝：全sourceと二つの全A軌道のmod-H affine kernel

2026-11-07。主義務G、producer plane_collision_review。
S⊂P4はunique hyperplane H内のquadric cone、π|S=0、一般rank四、実Dは既約被約quinticとする。
本稿は必要全sourceの有限還元であり、全Poisson支持又はG非優越はまだ証明しない。

## 1. 全六十lift、Euler五冗長性とabelian conormal

座標を `(x0,x1,x2,z,w)`、`H=(w=0)`, `q=x0x2−x1²` とする。
vertexはz方向で、I_S=(w,q)はprime complete intersection。
Koszul/Euler列から全projective sectionは

`Π=w L+q A`

で表せる。Lは全五十linear-coefficient lift、Aは全十constant skew係数。
E(-3)はacyclicなのでsectionを落としていない。
実150係数でraw source rank60、vertical rank24、合併rank79を再構成し、
projective section dimension55、Euler冗長性五を確認した。
冗長性はL→L+Euler∧constant vectorであり、Aはこの表示で一意である。

nonplane multiplicity二・Frobenius first-jetからgeneric conormal idealはabelian。
prime CIのI²はI-primaryなので、全三次bracket `{w,q}` はI²へ入る。
次数三のI²はw²linear及びwq constantの六係数で尽くされる。
全raw60とmembership六係数の35×66行列はrank17、完全kernel49。
membership係数のみを消すとraw60への全十一線形条件となる。
allowed Aの像は厳密に七次元で、`A^{w,xi}=0 (i=0,1,2)` だけを課す。
残りの条件はLの全sourceに保持した。

## 2. v=0又はC=0の枝はHを因子に持つ

`v=ΣA^{w,j}∂j=λ∂z`、`b_i=A^{xi,z}`、base3×3 skew部分をCと書く。
λ=0ならΠ|HはΛ²THに入り、H上でπ²=0、従ってH|D。
λ≠0でもC=0ならAは `∂z∧(−λ∂w−Σb_i∂xi)` とdecomposableであり、
`D|H=q² Euler_H∧A²=0`。やはりH|D。
これらは全L・全係数の枝であり、既約Dという仮定によってのみ除く。

## 3. 全Aの実normalizationと二軌道の完全cover

λ≠0ではnormal shear

`xi'=xi−(b_i/λ)w`, `z'=z`, `w'=w`

で全b_iが消える。逆は符号を戻すだけでdeterminant一。
qはw*linearだけ変わるので、その差は全Lへ吸収できる。
続くvertex scaling `z'=z/λ` で `A^{w,z}=1`。
これは未知Lを制約せず、λ=0の補集合は§2で保持済みである。
checkerは分母をλ²で払った全matrix恒等式とinverse shearを普遍多項式として検査した。

`C=C01∂0∧∂1+C02∂0∧∂2+C12∂1∧∂2` とすると、
SO(q)作用下のinvariantは `Δ=4C01*C12−C02²`。
dual linear formは `ell_C=C12*x0−C02*x1+C01*x2` で、非退化dual qに関するnormがΔである。
SO3(C)は非零normの各levelと非零isotropic coneにそれぞれ推移的。
非isotropic点では直交basis、isotropic点ではv,uをB(v,u)=1に取りuをnullに調整したWitt basisを作れば、
全点に実際の可逆orthogonal matrixが得られる。逆はQ⁻¹MᵀQである。
特定のGram pivotの零点を捨てるrational chartを仮定しない。

Δ≠0ではtensor scalar σをσ²Δ=−1に選び、変化したλをvertex scalingで再び一にする。
複素数体上でσは必ず存在する。
Δ=0,C≠0は非零nilconeの一軌道なので、残る全Aは次の二つで尽くされる。

```
Asemi=∂w∧∂z+∂0∧∂2,
Anil =∂w∧∂z+∂0∧∂1.
```

これは二つの任意標本による一般性の推論ではない。
Δ=0、C=0、λ=0の全補集合を以上で明示した。

## 4. 全mod-H Jacobiのaffine kernel

各Aを固定したまま、全五十Lを未定元として全350 projective Jacobiを構成した。
w=0への制限ではL²項が全て消え、全式はLについてaffine-linearとなる。
全conormal十一式も同時に課した。
両軌道ともaffine system rank20、全rawL affine dimension30を得た。
particular point、完全五十行kernel、全350のmod-w係数を保存した。

Eulerの五次元homogeneous方向が全kernelに入り、非零5minorに基づくliteral五座標gaugeを課すと、
各projective affine sourceは厳密に二十五変数となる。
全particular/kernel代入零とgaugeの完全性を確認した。
full Jacobiをこの二十五源で解いたのではない。これは全必要条件を保つsupersourceである。

## 5. 消費済みnormalization次元の復元

Sを保存するGL5 matrixは、baseのGO(q)、vertex線形変更、normal shear、w scalingのblock形を持つ。
次元は `4+3+1+4+1−1=12`（最後はglobal scalar）である。
AutS×tensor scalarは十三次元。
allowed AのΔ≠0 stratumは七次元、Δ=0,C≠0 stratumは六次元で、§3により各全軌道。
そのnormalized stabilizerはそれぞれ六・七次元となる。

今後、normalized Aにおける**full Poisson** fiberのdimensionをdと証明できれば、
全固定S affine sourceはd+7又はd+6以下。
AutS十二とprojective scalar一をfiberへ戻したPGL像上界は

`18+d`（Δ≠0）、`17+d`（Δ=0,C≠0）

となる。従ってG≤24に十分なのはd≤6又は7である。
今回得たmod-H supersource dimension25を、そのfull Poisson boundと取り違えない。

## 6. 有限保存物と次の構造課題

- `code/classification_quadric_cone_zero_source.py`
- `research/classification/jobs/quadric_cone_zero_source.json`
- `reports/CLASSIFICATION_QUADRIC_CONE_ZERO_SOURCE_20261107.json`

事前gate後、Fraction計算は1.757秒で終了。GB、jet反復、追加normal-form点列は用いていない。
全sourceと両affine kernelは不変保存する。
次の候補は二十五変数の幾何学的Poisson–Ore分解であり、full nonlinear支持を
無根拠な長時間GBへ渡すことはしない。Gはこのquadric-cone枝についてなお未完了である。
