# 残h5実軌道：一つのforced canonical ambient契約と純構造

2026-11-08、PRIMARY `E_multiple_hyperplane`。**契約準備のみ、未GO・未実行。**
full h5 sourceの既知閉和＋一軌道への全還元を入力とする。

## 同じ一点を使う理由

`UNIMODULAR_LOG_POLE_WHOLE_BETA_ZERO_AND_SINGLE_ORBIT.md` §5の全source正逆により、
残sourceは全translation/最高GL4/tensor/unitを戻して一つのactual orbitである。
そのforced pointは順序(xy,xz,xw,yz,yw,zw)で

```
P_*=(x,−x²,xz−1,xy−1,−yz+w−y,−xw+xy+y+z).
```

旧boundary reportはこの全sourceの24/6低次係数kernelと全残式を保存するが、
このforced pointの全350×150 ambientD／全51group方向の証明書は含まない。
同一点の整数liftはすでにscale一であり、全sourceから再入力出来る。
既存ambient engineの変更や別sourceのoperator再生成は不要。

## 純構造

最高QはX∧Yで

```
X=x∂x−y∂y+(y−w)∂w, Y=−x∂z+z∂w,
[X,Y]=Y, trX=−1,trY=0.
```

全P_*にはaffine Poisson場

```
Z=x∂z−z∂w−∂y
```

が有り、全六Lie成分は零。
standard volumeによるcomplementary二形式ω（matrix inverseとは符号が異なる）を使うと

```
i_Zω=y dx+x dy+dz=d(xy+z).
```

よって第一添字sharp規約ではZ=−P_*#dHであり、quadratic noncritical potential −H=−xy−zのHamiltonianでもある。
homogeneous Zは(h,x,y,z,w)上でpartition(3,2)、chains h→y、x→z→w。
これはactual parentを同定した証明ではない。
特にHがquadraticであることをlinear-Hamiltonian最高部定理へ取り違えない。

## ONE full ambient判定の最弱な矢印

同じP_*一つだけで、全150quadratic lift、全350projective Jacobi、全25Euler列、
全25GL5列とtensor scalar列を作る。全strict Jacobi／Pf／divも元sourceから戻す。
fullD、EulerV、fullMの各完全kernel/freeidentityと非零rank minor、DM=0を保存する。

r=rankD、m=rankMとする。Euler rank24とscalarを除いて

```
projective tangent dimension=125−r,
actual orbit dimension=m−25.
```

r+m=150なら実軌道の次元下界と全接上界が一致し、P_*はsmooth、
orbit閉包が唯一の全空間成分となる。新規性／既成分との同定はその後の独立監査へ分ける。
r+m<150なら、同一点の局所経路は正直に非閉鎖として止め、別point／higher jet／minor列を追加しない。
実既知親への幾何的曲線又は全係数同定へ戻る。

この契約は一つのwhole sourceのforced orbitの判定であり、generic性を期待して選んだ点の列ではない。
純actual parent curveが先に得られればambient自体を実行しない。
resource capは全job180秒／3GiB、失敗時も元source／全完了データとfailureをimmutable保存する。
