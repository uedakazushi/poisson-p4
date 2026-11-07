# quadric cone上の全rank零枝：Poisson–Ore分解とPGL像二十三上界

2026-11-07。主義務G、producer plane_collision_review。
**命題。** quadric cone S⊂P4上で零となり、実Dが既約被約quinticである全rank四Poisson構造の
PGL5飽和像のdimensionは二十三以下である。
root/別作者の独立受理前はproducer証明として扱う。
clean、node又はvertical成分に関する仮定は置かない。全E帰属やrank四不存在は主張しない。

## 1. 全sourceと正規化の入口

`QUADRIC_CONE_ZERO_SOURCE.md` の全Koszul/Euler source、abelian conormal必要条件、
全A二軌道coverを入力とする。
座標は `(u0,u1,u2,z,w)`、q=u0u2−u1²、S=(w=q=0)。
λ=A^{w,z}=0又はbase skew C=0ではH=(w=0)がDを割るので、既約Dの仮定から除く。
残る全Aは

`A=∂w∧∂z+C`、`C=∂0∧∂2` 又は `∂0∧∂1`

へ実AutS×tensor scalarで移る。
allowed A stratumのaffine dimensionは七又は六、AutSは十二次元である。
いずれも二つの任意点から一般性を推測したのではなく、Δ≠0及びΔ=0,C≠0の全軌道である。

## 2. 全二十五変数の可逆Ore bridge

X=C*dq、u=(u0,u1,u2)とする。
全mod-w Jacobi/conormal affine kernelをEuler五冗長性で割ると、次の全二十五変数表示を取れる。

```
L_base=w D,                Dは任意constant skew3×3,
L^{w,base}=−X+w a,         a∈C³,
L^{w,z}=b·u+c z,           b∈C³,c∈C,
L^{z,base}=M u+z d+w e,    M任意3×3,d,e∈C³,
Π=w L+q A.
```

この表示の25parameter方向と全Euler5方向が、保存済み全30raw kernelそのものを張ることを確認した。
particular pointの差も同じkernelに入り、30×50のexplicit recovery matrixで全係数を逆に回収した。
従って全kernelの一部分源を選んだのではない。

chart w=1でprojective bivectorを正確にdehomogenizeすると

```
P=B(u)+(Y(u)+z T(u))∧∂z,
B=q C+u∧X+D−u∧a,
T=a−d−c u−X,
Y=−u q−u(b·u)−M u−e.
```

translation変数t=a−dへの変更は可逆であり、Yの十五係数にも追加制限は無い。
全projective350係数のJacobi条件は、次の全条件と両方向に同値である。

`Jac(B)=0`, `Lie_T B=0`, `Lie_Y B=Y∧T`。             (1)

checkerはglobal E∧Jac(Π)とaffine Jacobiのidentity、さらに(1)へのidentityを全25変数上で直接検証した。
原点近傍だけのjet同値ではない。

## 3. BとTの全支持は三次元以下

Jac(B)の空間次数二・一の係数はD,aの六係数に対するrank三の線形式である。
それらは正確に `a∈image C`, `D=δ C` を与える。
この三次元線形空間へ全Jac(B)を戻すと恒等零となるので、Bの全Poisson支持である。

### Δ≠0軌道

`a=(p0,0,p1)` と置く。
T=t−c u−XがBのPoisson derivationとなる全条件は

```
t=((c+1)p0,0,(1−c)p1),
c*(δ−p0*p1)=0.
```

c=0ではp0,p1,δの三変数、c≠0ではδ=p0p1でp0,p1,cの三変数である。
特に固定した非零cではbase dimensionは二以下。

### Δ=0,C≠0軌道

`a=(p0,p1,0)` と置く。
全条件は

```
t1=p1−c*p0/2, t2=0, t0自由,
c*p1=0,
−4cδ+c*p0²−2p0p1+2p1t0=0.
```

c≠0ではp1=0、δ=p0²/4で、p0,t0,cの三変数。
c=0では `p1*(t0−p0)=0`、δ自由なので、各閉支持は三次元である。
固定した非零cでは二次元以下。

全条件は定数pivotによるtの線形消去と、残った全係数のliteral代入で得た。
パラメータを反転してc=0、p1=0、nilpotent又は零baseを捨てていない。
両軌道の全(B,T)支持は、したがってdimension三以下である。

## 4. 全Y fiberの三角形rank上界

Yの未知十五係数を、radial quadraticのb三係数、linear M九係数、constant e三係数に分ける。
fixed cubic項は−u qである。
最後の(1)はYについてaffine-linearで、空間次数三、二、一の順に解ける。
各段階の未知係数に掛かるdiagonal matrixは、B0=qC+u∧X、T1=−X−c uだけで決まり、
低次数のa,δ,tには依存しない。

従って各段階の全nullityの和は、**任意の全base値**におけるY affine fiberのdimension上界である。
低次数coupling又はfixed cubic項は右辺を変え、fiberを空又は小さくすることはあるが、
この逐次nullity上界を大きくしない。最高次数又は最後の定数方程式を追加しても上界は保たれる。

全diagonal matricesをQQ(c)上で計算し、generic非零minorと完全kernelを保存した。
minorの全因子をexactに分解したところ全根は有理数であり、全例外をそのまま再評価した。
generic RREFの分母に依存せず、literalなminor非零性をrank下界として使うので、
以下以外のcを暗黙に捨てていない。

### Δ≠0軌道

| block | 大きさ | generic rank | 一つのgeneric minor |
| --- | --- | ---: | --- |
| b | 30×3 | 2 | −1 |
| M | 18×9 | 8 | `(c−2)(c−1)²(c+1)³(c+2)²` |
| e | 9×3 | 3 | `(c−1)(c+2)²` |

generic Y fiberは `1+1+0=2` 以下。
全例外c=−2,−1,1,2では三rankが(2,7,3)、fiberは三以下。
constant blockの表示minorが零となる場合も、別の非零minorでrank三を確認した。
c=0では(2,8,3)のまま、fiber二以下である。

generic base三次元＋fiber二、c=0のbase三＋fiber二、
非零例外cのbase二＋fiber三の全てで、normalized full sourceはdimension五以下。

### Δ=0,C≠0軌道

generic ranksは同じく(2,8,3)、minorはそれぞれ1,c⁸,c³。
全例外はc=0のみで、rankは(2,7,3)、fiber三以下。
generic normalized sourceは三＋二で五以下、c=0は三＋三で六以下となる。

## 5. 全PGL像への復元

normalized full Poisson fiberの上界は、Δ≠0でd≤5、Δ=0,C≠0でd≤6。
§1のA軌道を戻した全固定S affine section源は、各々d+7又はd+6で、ともに十二以下。
AutS十二次元とtensor scalar一次元は、PGL5×非零sourceのfiberに十三次元の自由な部分を与える。
個別tensorのstabilizerが零とは仮定しない。
従って両軌道の全PGL像は

`24+12−12−1=23`

以下である。λ=0又はC=0の補集合は既約Dでは存在しないので、これで全quadric-cone rank零枝を覆う。
未知二十五〜二十七次元G成分の一般元から排除される。

## 6. 不変証明書と射程

- 入口: `code/classification_quadric_cone_zero_source.py`、同名job、
  `CLASSIFICATION_QUADRIC_CONE_ZERO_SOURCE_20261107.json`、`QUADRIC_CONE_ZERO_SOURCE.md`。
- Ore全有限計算: `code/classification_quadric_cone_zero_ore.py`、
  `jobs/quadric_cone_zero_ore.json`、`CLASSIFICATION_QUADRIC_CONE_ZERO_ORE_20261107.json`。

Ore計算は事前gate後1.502秒、全項目PASS。25-variable GB、jet反復、任意点probeは無い。
二度のpreflightはSymPyのmatrix-assignment/symbol-list APIの補正で、報告を作る前に停止した。
profile又は数学の範囲を変えて再計算していない。
生成器をimportしない独立replayと幾何的coverの別読解後に限ってglobal主張へ昇格する。
E帰属、新成分数、rank四不存在、他の曲面又はG全体の完成は結論しない。
