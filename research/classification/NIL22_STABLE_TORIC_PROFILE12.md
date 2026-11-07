# stable_three 全磁気源の Cartan 正常形と一回の構造的 toric confluence

2026-11-07。E_multiple_hyperplane。
全六磁気源を三変数の完全源へ送る実変換を構成した。
その後、root が許可した一つの (1,2) profile を全 linear r 上で調べたが、
一般 r では不整合である。残る r の真の閉補集合を明示し、toric route を停止する。
T01 非所属も E 完了も主張しない。

## 1. 全六磁気からの実 Cartan bridge

元座標を `(h,q1,q2,x,y)` とし、次の固定 PGL5 変換を行う。

```
H=3h+x,
Z=−q2/3,
W=(q2−q1)/3,
U=3h−y,
V=3h−x+y.
```

H 行を Euler-vertical 補正で消し、tensor を9倍すると、affine H=1 上で
全六磁気を保つ源は

```
ell=1+U+V,
{Z,U}=ell*U, {W,V}=ell*V,
{Z,V}={W,U}={U,V}=0,
{Z,W}=−V*Z+U*W+R(U,V)
```

となる。R は全二次式で、元の六磁気係数と可逆な線形変換で結ばれる。
R の U,V,UV 係数をそれぞれ cU,cV,cUV とすれば、

`Z'=Z−cV−cUV*U`, `W'=W+cU`

でその三項が消える。全パラメータで可逆であり、磁気係数は反転しない。
完全な残余源は

`R=a0+a1*U²+a2*V²`

である。三 a は元 m に関する独立な線形形式であり、係数行列の rank は3。
その全表示と全六括弧の恒等式は V2 報告に保存した。
従ってこの三変数源への全支持認識は、元の全六磁気源への認識を与える。

## 2. 最初の (1,1) profile の rank4 障害

H,U,V を固定した toric 座標を
`S=ell+tZ`, `T=ell+tW` とする。
toric period `kij(t)` は t に関し regular である。
実際、target 側で `{Li,Lj}/(Li Lj)` は空間変数によらない period に等しく、
全 Li の極限が非零である点へ評価すれば regular 性が従う。

上の対角 anchor を極限に得るには、period の一次係数が

```
kUS,1=−1, kUT,1=0,
kVS,1=0,  kVT,1=−1, kUV,1=0
```

となる必要がある。一方 `{Z,W}` の t⁻1 係数の消滅は

`kST,1*ell−U+V=0`

を要求する。ell の定数項から kST,1=0 となり矛盾する。
一般の anchor 係数で書けば同じ比較は二列の比例を強制するので、rank4 を作れない。
磁気項をどの値へ置いても、この障害は変わらない。

## 3. 一回だけ調べた全 linear-r (1,2) profile

次の全 linear form を保持した。

```
r=r0*H+r1*U+r2*V,
L3=ell+t*Z,
L4=ell+t*Z+t*r+t²*W.
```

toric 座標は `(H,U,V,L3,L4)`、変換行列式は `t³`。
逆行列の pole は高々2、bivector の pole は高々4なので、regular period の
Taylor 係数0,...,4だけで全負冪と端点が完全に決まる。
六 period の全30係数を残した。これは Jacobi jet を有限次数だけ合わせる方法ではなく、
t≠0では常に実 toric Poisson tensor である族の、厳密な有限端点問題である。
成功した場合には period を次数4多項式に取れば全 t の実 Poisson 曲線になる。

全負冪消滅と全六端点括弧を `QQ(r0,r1,r2)` 上で一度に比較した。
30 period 係数、三磁気変数、右辺を合わせた全行列は81×34、period rank21、
augmented rank22 であり、compatibility は `−1=0` となった。
磁気の特殊値を選んでも、この generic-r 不整合は解消しない。

## 4. 不変な証明書と、捨てていない r 補集合

V2 は元の全81行、空間 monomial/次数 labels、全非零 RREF 行、
元の22行による exact functional を保存する。
functional の唯一の分母は `r0−r2` であり、元の全34列に対して

`sum(weight_i*row_i)=(0,...,0,1)`

を直接確認した。分母を払った同じ多項式恒等式から、任意のこの profile の解には

`r0−r2=0`

が必要と分かる。従ってこの真の閉補集合は未分類のまま保持する。
四番目の hyperplane の rescaling による `r→r+c*ell` を考慮すると、
そこでは r は U の倍数に相当するが、その枝の実行は今回行っていない。
一般 r の不整合を全 r の非存在へ拡張しない。

指定された toric route はここで停止した。(1,3)、高い period 次数、
特殊 r の scan、別 frame farm は行っていない。
別の実変形又は別成分への帰属は引き続き未解決である。

## 5. 保存物

- `code/classification_nil22_stable_toric_profile12.py`
- `research/classification/jobs/nil22_stable_toric_profile12.json`
- `reports/CLASSIFICATION_NIL22_STABLE_TORIC_PROFILE12_20261107.json`（最初の不変実行、1.025秒）
- `reports/CLASSIFICATION_NIL22_STABLE_TORIC_PROFILE12_V2_20261107.json`（完全な線形証明書、1.265秒）

V2 追加は行列・RREF・functional の保存だけで、profile、次数、全源の範囲を変えていない。
各実行前に job gate を通した。
これらは producer 証明書であり、生成器を import しない別作者 replay は未実施。
