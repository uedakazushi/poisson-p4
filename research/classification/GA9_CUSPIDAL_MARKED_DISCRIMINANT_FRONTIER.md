# Ga9平面quartic24成分：残る二直線・三次三配置を含まない

2026-11-07、PRIMARY `E_reduced_nonnormal`。新CAS、点、法線又はjetを用いない手証明。
既知成分の閉退化因子に対する必要条件を全marked flagで閉じ、全源の親を一つ除外する。
他のambient親又は新成分性は主張しない。

## 1. 使う既知成分と全源

`GA9_CUSPIDAL_QUARTIC_COMPONENT.md` とその独立
`GA9_CUSPIDAL_QUARTIC_SEPARATE_REVIEW.md` で受理された
`R4_GA9_CUSPIDAL_PLANE_QUARTIC_24` をCとする。
全六磁気係数を持つ実sourceのPGL5飽和閉包がCそのものであり、実退化式は非零scalarを除き

```
D=h F_a, F_a=(hW+6x²)²−144x³(h+ax).
```

全磁気係数はDに寄与しない。任意のa、特に一般aについて、
H={h=0}とquarticは一点p=[h:x:W]=[0:0:1]で全交点重複度4を持つ。
一般aでのquarticのpにおける二次coneはh²、すなわちHの二重線である。
以下の全閉条件は一般aで成立すれば十分であり、source chartの分母やquartic既約性の境界も
properなflag像の閉性によりCの全体へ含める。

## 2. exact essential3の商平面とproper flag

実五次Dの方向微分写像のrankが正確に3である開では、そのkernel Kは二次元bundleである。
従ってDの本質三変数商平面P(V/K)は正則に復元でき、PGL5の座標混合によって
この商平面のmarked幾何が失われることはない。
一般sourceからrank3のtargetへ近づく際も、このrank3開内でkernelと商を同一視できる。

H、quartic F、p∈Hのflagを伴わせる。D=HFという射影乗法写像、
H上の点pのincidenceはproperである。
target Dが相異なる二単純超平面H1,H2と既約三次からなる場合、
極限の親Hは正確にH1又はH2のどちらかである。
既約三次の一部分へHが紛れ込むことや、別の超平面が生じることはない。
各選択に対しF=D/Hは一意である。

親sourceではF|Hが四乗linear formなので、この条件もproper乗法像で閉である。
targetのF|Hは零ではない（DはHを重複して含まない）。
従ってこれが四乗でなければ、そのH選択は排除される。
四乗ならpは一意であり、次節の条件はこのpに適用される。

## 3. marked binary discriminantの三重根条件

flagをH={h=0},p=[0:0:1]へ移す。
F|Hがx⁴の非零倍、かつpの二次coneがh²の倍であることから、Fは必ず

```
F=a0 h²W²+hW B2(h,x)+C4(h,x)
```

と書ける。B2はbinary quadratic、C4はbinary quarticである。
この表示はa0=0の場合も含み、a0を反転していない。
二次coneがh²の倍である条件はlinear係数条件なので、cone全体が零となる境界も含め閉である。

binary quartic

```
G4=B2²−4a0 C4
```

を取る。flag座標の取り替えはbaseの可逆linear変換と
`W→dW+L(h,x)`（d≠0）で生成される。
後者をFへ代入すると

```
a0_new=a0 d²,
B2_new=d(B2+2a0 hL),
C4_new=C4+hL B2+a0 h²L²,
G4_new=d² G4.
```

base変換でもG4は非零scalarとbinary linear変換を受けるだけである。
従って「G4が三重根を持つ、又はG4=0」はflagに内在的である。
これはl³mの射影乗法写像P1×P1→P(Sym4 C2)のproper像のaffine coneであり閉条件である。
零quarticを除外していない。

親の全F_aでは厳密に

```
a0=1, B2=12x²,
C4=−144hx³+(36−144a)x⁴,
G4=576x³(h+ax).
```

従って全一般sourceが三重根条件を満たす。
§2のproper flagと合わせ、この必要条件はCの全boundaryのうち
exact essential3でDが二単純H+既約三次となる全点にも適用される。
quadratic leading coefficientの零点、三重根と第四根の衝突、H及びpの座標chartの境界を
いずれもこの閉性の中に保持する。

## 4. node_flex_tangent_through_flex

三次は `C_N=Y²Z−X²(X+Z)`、二直線はZ,X+Zである。
親HをZと選ぶとF=(X+Z)C_N、F|Z=−X⁴なのでp=[0:1:0]。
Y=1のlocal座標X,Zで二次coneは

```
q2=Z(X+Z).
```

これはH²=Z²の倍でないため、§3の閉二次cone条件を破る。
もう一つのH=X+Zを選ぶとF=Z C_Nで、そのHへの制限はZ²Y²である。
異なる二根を持つので四乗linear formでなく、§2で排除される。
従ってこの全磁気sourceにCの点はない。

## 5. cusp_flex_secant_through_flex

三次は `C_K=Y²Z−X³`、二直線はZ,X−Zである。
H=ZならF=(X−Z)C_K、制限は−X⁴、同じp=[0:1:0]で

```
q2=Z(X−Z)
```

となり、再びH²の倍でない。
H=X−ZならF=Z C_Kの制限はZ²(Y−Z)(Y+Z)で、四乗でない。
この全磁気sourceもCの外にある。

## 6. cusp_cusp_tangent_cusp_secant

三次は同じC_K、二直線はY,Y−Xである。
H=Yを選びh=Y,x=Y−X,W=Zと書くと

```
F=x[h²W−(h−x)³]=h²xW−x(h−x)³.
```

F|h=x⁴なので唯一のp=[0:0:1]は許され、二次coneは零である。
しかしa0=0、B2=hx、C4=−x(h−x)³から

```
G4=h²x².
```

これは相異なる二つの二重根を持ち、三重根を持たない非零quarticである。
従ってa0=0の境界を含む§3の閉三重根条件を破る。
他のH=Y−X（x=0）を選ぶとF=Y C_Kの制限はh³(W−h)で、四乗でない。
よってこの全磁気sourceも、I=0等の磁気境界を含めCの外にある。

## 7. 結論と残る親

三対象のDと全marked配置は磁気六係数によらない。
上の除外は特定の磁気一点、固定したsource chart、又はtangent upperに依存しない。
全PGL飽和源に同じ結論が成り立つ。

従って既知Ga9 cuspidal plane-quartic24成分は、残る平面1+1+3三配置のいずれの親でもない。
これはT36/T37/T38等の親の除外、新ambient成分の構成、又は全source被覆ではない。
独立読解の対象はproper flag/essential quotient/全二H選択とbinary三重根の閉性である。
