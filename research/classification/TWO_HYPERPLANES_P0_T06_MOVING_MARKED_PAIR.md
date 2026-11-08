# p=0 全五源と T06 cubic の実 moving marked-pair class

2026-11-08、PRIMARY `E_reduced_nonnormal`。旧四法線・高次 jet を使わず、
全源 normalizer と同じ cubic の unipotent automorphism から一つの actual MPS 親 class を導く。
以下の有限 matching は未実行で、陽性の帰属を先取りしない。

## 1. 元の全五源の厳密な一 parameter normalizer

`TWO_HYPERPLANES_CUBIC_SUPPORT.md` §5.2 の p=0 必要十分五源を使う。
自由変数は `n2,n3,n5,n8,n10` で、`S=h+x` とする。

```
A=n10*n3+n5+n10²/2, B=n2+n3*n8,
U=u+A*h-B*x, V=v, W=w-n10*h-n8*x.                            (1)
```

これは det=1 の全三角 PGL、逆は `u=U-A*h+B*x,v=V,w=W+n10*h+n8*x`。
source parameter の分母を持たず、exact cofactor `V0=v∂u+S∂v+h∂w` を保つ。
変更後は `l=n3` として

```
{x,u}=xv, {x,v}=xS, {x,w}=hx,
{u,v}=hu-v²/2+l*(-hv+Sw),
{u,w}=hu-w²/2+l*h*(w-v),
{v,w}=h*(v-w)-l*hx.                                         (2)
```

全 l で実三次は

```
F0=[h*(v-w)²+x*w²-2hxu]/2.                                 (3)
```

u の係数 hx と残部は互いに素なので既約。
全5 source は (2) の実 PGL 飽和閉包へ入る。`diag(1,1,t²,t,t)` pushforward と tensor t^-1 は
l を tl へ送るので、l 非零は一つの projective orbit。これは既 λ=1 点を別の点へ置き換える操作ではない。
`Z0=(2v+w)∂u+(3h+2x)∂v+3h∂w` は (2) の全 l で同じ実 affine modular field で、
`Pi_l=Pi_0-l V0∧Z0`。ただし以下はこの twist が他の component に留まると仮定しない。

## 2. 真の T06 source と同じ cubic

受理済み raw T06 の全四 parameter `(p0,p1,t0,t1)` source をそのまま使う。
source tensor は四 parameter に線形で、全値で元の全350 projective Jacobi を満たす。
`L=2p0+p1` として実五次の三次因子は

```
L*x0*x1*x2+t0*x2*x4²-t1*x1*x3².                             (4)
```

`t0=L,t1=-L` は C=L*(UVW+UX²+VH²) を与える。
ここで `(x0,x1,x2,x3,x4)=(W,V,U,H,X)`。
`(U,V,X,H,W)=(h,x,v-w,w,-2u)` なら C=2L*F0。
ただし marked hyperplanes は source の H,X であり、target の h,x と違う。

この cubic の actual unipotent automorphism は

```
X->X+t*V, H->H+s*U,
W->W-2t*X-t²*V-2s*H-s²*U.
```

t=s=1/e を取ると source へ送る actual frame は

```
x0=-2u-2v/e-(h+x)/e²,
x1=x, x2=h,
x3=w+h/e, x4=v-w+x/e.                                     (5)
```

det は非零定数、逆は

```
h=x2, x=x1, w=x3-x2/e,
v=x3+x4-(x1+x2)/e,
u=-x0/2-(x3+x4)/e+(x1+x2)/(2e²).
```

従って正規化した source H,X は e=0 で h,x へ移り、同じ cubic を保つ。
これは旧 fixed flag の再試験でなく、intrinsic marked pair を動かす一つの actual PGL class である。

## 3. 先験的に有限の全 period/edge matching

(5) と逆は e の pole 高々二。二次 tensor の pushforward は
座標代入の pole 高々四と wedge Jacobian の pole 高々四から、一様に pole 高々八を持つ。
Euler 商は `C(P)=P-E∧div(P)/5` の係数線形射影で取り、同じ bound が成り立つ。
元の source 四係数の Euler 商への写像は定数の rank 四であり、定数 QQ 左逆を持つ。
従って (5) で動かす class の regular tensor 曲線から source parameter を戻すと、
その四 Laurent parameter の pole は高々八。
逆方向の pole bound も八なので、parameter の e^k, k>8 は endpoint に寄与しない。

以上から、この derived fixed-frame class における全 meromorphic regular limit は
各四 parameter の `e^-8,...,e^8` の17係数、計68未知数だけで厳密に判定できる。
source は四係数に線形であり、parameter 積の独立緩和をしていない。
全負冪係数零と e^0 全150 Euler 商が (2) に一致する条件は、一つの QQ 線形系である。
(2) は l に線形なので、二本の RHS を同時に保持すれば全 l を判定できる。

陽性なら保存する実 parameter profile を元の真 T06 source へ直接代入し、
全150逆・負冪・全350・Pf・全 l endpoint を検証する。閉性と (1) の全逆で全五源が覆われる。
陰性なら完全余核 functional から l の正確な有限 exceptional locus を保存し、
この一 moving-frame class の陰性までに留める。全 T06 不存在や全 known-parent 陰性へ上げない。
別 profile の増殖、Taylor bound 増加、別点又は高次 normal 計算は行わない。

## 4. 既データが与える別の次元制限

既存 `CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_NONZERO_MODULUS_20261107.json`
の同じ l=1 点は全 D102・実 orbit/source47、Euler24。
従って ambient projective tangent 上界は23、全24/25成分はこの点を含めない。
ただし同次元以下の既知 source は別の actual 認識/排除を要し、
これだけで新成分又は全 catalogue 外としない。上の matching は dim23 の T06 を一つ扱う。
