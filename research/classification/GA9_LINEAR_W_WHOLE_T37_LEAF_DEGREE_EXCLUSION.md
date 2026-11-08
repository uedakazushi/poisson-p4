# Ga9 linear-w は whole T37 の外：四次成分の全 leaf は次数九

2026-11-08、PRIMARY E、純手証明。source／specialization は `T37_WHOLE_CLOSED_LOW_DEGREE_LEAF_SURFACE_THEOREM.md`。新しい operator、point、rank、minor、normal/profile 又は jet は作らない。

## 1. 対象の全 unit 軌道

既受理 `GA9_LINEAR_W_FIXED_DIVISOR_MAGNETIC_NORMAL_FORM.md` により、退化因子

```
D=−4hF4, F4=4h³w−8h²y²+y⁴
```

を持ち、唯一の singular plane S=(h=y=0) 上で一般 rank二の全 source は、実 fiber GL2／全六 magnetic shear とその全逆で次の一つの actual23 軌道へ戻る：

```
Pxy=4h²−y², Pxz=−hx+2yz+4w²,
Pxw=16hy−4wy, Pyz=hy, Pyw=0, Pzw=−4hw+4y²,
h row=0.
```

磁気零の別軌道へ非所属を外挿しない。F4 は primitive w-linear で絶対既約。

## 2. 四次 component 上の全 affine leaf

h=1 では F4=0 が w=2y²−y⁴/4 を与え、三変数 (x,y,z) 上の括弧は

```
{x,y}=4−y², {z,y}=−y,
{x,z}=−x+2yz+4(2y²−y⁴/4)².
```

K=yx+(4−y²)z と置けば

```
{x,K}=4(4−y²)w², {z,K}=−4y w².
```

従って

```
C=K−16y⁵/5+4y⁷/7−y⁹/36                                (1)
```

は Casimir。微分された最後の polynomial は 4w²=16y⁴−4y⁶+y⁸/4 に一致する。y と 4−y² は同時に零でないので、π の rank はこの全 affine open で厳密に二。

全有限 κ に対し C=κ は integral A² である。実際

```
T=−x/4+yz/4,  K=yx+(4−y²)z,
x=(y/4)K−(4−y²)T, z=K/4+yT
```

は C[y] 上の determinant 一の可逆変換で、K=C+16y⁵/5−4y⁷/7+y⁹/36。よって (y,T,C) は affine A³ の全 polynomial coordinates である。各二次元 leaf は (1) の fibre、全 fibre は連結である。

## 3. Leaf closure の次数を全て九と計算する

F4 の birational resolution を

```
Y=P_P1(O(4)⊕O⊕O),  L=tautological class, f=P²-fibre class
```

で与える。base coordinates を [s:t]、fibre coordinates を (α,β,γ) とし、α の class は L−4f、β,γ の class は L。actual map は

```
(h,y,w)=α*(s⁴,s³t,2s²t²−t⁴/4), x=β,z=γ.                 (2)
```

general base point で rational quartic curve の normalization parameter は y/h=t/s だから、(2) は h≠0 上で birational。exceptional divisor α=0 の contraction 又は cusp fibre が general leaf の degree を増減させると仮定しない。

(1) の fibre C=κ の closure は、(2) へ戻して αs⁹ を払った literal equation

```
s⁴t β+s³(4s²−t²)γ
 −α[(16/5)s⁴t⁵−(4/7)s²t⁷+t⁹/36+κs⁹]=0.                (3)
```

である。全項は class L+5f。β,γ 係数の共通因子は s³、しかし α 係数は s=0 で −t⁹/36≠0。従って三つの fibre-linear coefficients に共通 base zero はなく、(3) は全有限 κ で primitive な相対 hyperplane、integral divisor。α=0 又は s=0 の二次元 fixed component を含まない。

Y で L³=4,L²f=1。像への map は (3) の一般点で birational なので、leaf closure の actual projective degree は

```
L²*(L+5f)=4+5=9.                                         (4)
```

これは二つの方程式の Bézout 数をそのまま degree とした計算ではない。全 exceptional／infinity 部分を含む一つの integral divisor (3) の push degree であり、全有限 fibre に同じ値が成り立つ。

## 4. Whole T37 の必要条件との矛盾

もし canonical Ga9 が whole closed T37 に属せば、共通 source theorem の DVR/Hilbert 議論により、F4 component の一般 smooth rank二 point p は degree≤5 の invariant surface Z に属する。

p は h≠0 と取れる。そこで π の全 leaf は integral C=κ fibre。irreducible invariant surface Z は p の二次元 analytic leaf を含むので、その Zariski closure (3) を component に持つ。同じ dimension 二で integral な Z 自身がその closure となり、(4) により degree九。degree≤5 に矛盾する。

従って Ga9 unit 全23軌道は whole T37 の外。source の all zero magnetic boundary の非所属、新 componenthood／smoothness／成分数はこのノート単独では主張しない。既知他親の除外と全候補被覆の統合は別の義務として保つ。
