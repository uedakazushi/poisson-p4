# T21の本質三変数・二直線三次境界は二節点接線型だけ

2026-11-08、PRIMARY `E_reduced_nonnormal`。純証明であり、新CAS・anchor・jetは無い。

## 全源入口とT21の閉必要条件

実Dが相異なる二単純超平面と既約plane cubicからなり、本質変数数が正確に3とする。
`TWO_LINES_PLANE_CUBIC_SUPPORT.md` の全源還元を使い、baseを(h,x,y)、fiberを(z,w)、
二直線をh,xとする。全Poisson代表は

```
P=U∧∂z+V∧∂w+(−az−bw+R)∂z∧∂w,
det(U,V)=unit*x*C(h,x,y), U_h=V_h=0,
U_x,V_x はxで割れる。
```

Rは任意のbase二次式であり、U,V,a,bは全源に対するaccepted Saito還元から来る。
本稿は先にRを零とした切片へ主張を狭めない。

`TWO_LINES_CUSP_TANGENT_SECANT_GLOBAL_PARENT_EXCLUSION.md` §§2–3とその独立reviewにより、
T21全成分のこの種類のboundaryでは、Z_hをh=1代表の発散、X=P#dlog(x/h)として

```
Z_h+X は x=0 上でprojectively零、
Z_h+4X は h=0 上でprojectively零
```

が必要である。二条件はh,xの交換で交換され、proper factor flagで全境界に保たれる。
これはT21の一般divisor配置だけから推測したものではなく、全period/edgeを残す
完成されたT21 tensorの恒等式と、正確な二単純H因子flagの閉性から来る。

## Liouville恒等式から二つの三重接触を強制する

Pfaffianは非零定数を除きxCである。symplectic openでのLiouville volume
`dx dy dz dw/(xC)` がPに保たれることから、発散の規約を
`Z^j=Σ_i ∂_i P^{ij}` として厳密に

```
Z_h=−P#dlog(xC),
Z_h+X=−P#dlog C,
Z_h+4X=P#dlog(x³/C).
```

これらはhomogeneous rational fieldの恒等式なので、元の非零開から全係数に延びる。
全三場はfiber方向だけを持つ。したがってh=0又はx=0上でprojectively零なら、
そのhyperplaneの一般点では文字通り零である。

x=0上ではU_x=V_x=0である。
U_yとV_yが両方xで割れればdet(U,V)はx²で割れ、xCの単純因子xに矛盾する。
従って少なくとも一方はx=0の一般点で非零である。
第一の必要条件はその点で `C_y|x=0=0` を強制するから

```
C|x=0=c_h h³, c_h≠0.
```

Cは既約でxを因子に持たないので、c_hは零でない。
この議論はC≠0の一般点を使って必要な多項式恒等式を得ただけであり、
C=0上の未知Poisson点を排除したものではない。

h=0上ではPf=xCが一般点で非零だから、x,y,z,wの四次元行列Pは可逆である。
第二条件と上の恒等式から、h=0におけるx³/Cのx,y微分が零となり、
斉次次数を比較して

```
C|h=0=c_x x³, c_x≠0
```

を得る。hで割る操作もCの特殊点を捨てる操作もしていない。

## 唯一のmarked cubic型と残三配置への適用

二つの制限からCは

```
C=c_h h³+c_x x³+hx(d_h h+d_x x+d_y y)
```

となる。d_y=0ならbinary cubicは複素数体上で可約なので、d_y≠0。
従ってp=[h:x:y]=[0:0:1]の二次coneはd_y hxであり、pはnodeで、
h,xはその二本の相異なる節点接線である。
これはまさに `node_two_node_tangents` 型であり、その**全六磁気源**は
`PLANE_CUBIC_TORIC_COLLISION.md` と独立reviewの実族により既にMPSへ被覆されている。
逆に任意のmagnetic tupleがT21そのものに入るとは、この必要条件だけから主張しない。

残る三源をaccepted座標で書くと、

```
node_flex_tangent_through_flex:
    C=h y²−x(x−h)², C|x=0=h y²,
cusp_flex_secant_through_flex:
    C=h y²−(x+h)³, C|x=0=h y²−h³,
cusp_cusp_tangent_cusp_secant:
    C=h²y−(h−x)³, C|x=0=h²y−h³.
```

三つともc_h h³ではない。よって**残る三配置の全六磁気源はT21の外**である。
条件は全RとPGL飽和に一様で、磁気零も含む。
T36/T37/T38、別のambient成分、又は全E被覆への結論は追加しない。
