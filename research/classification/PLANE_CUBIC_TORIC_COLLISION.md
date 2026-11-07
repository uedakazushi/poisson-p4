# 節点三次と二本の節点接線：全磁気源の MPS 被覆

2026-11-07。主課題 `E_reduced_nonnormal`。
`TWO_LINES_PLANE_CUBIC_SUPPORT.md` §6 の未決四配置のうち
`node_two_node_tangents` の**全六パラメータ源**を、実際の多項式 Poisson 族により
既知 MPS40 成分の有限閉和へ含める。新成分数の追加はない。
他の三配置又は nil31 short r=2/3 の被覆はこの証明には含まれない。

## 1. SNC 点を含む補助の既約源

五つの独立座標を `(H,A,B,C,W)` とする。H-Casimir 代表において
任意の三変数斉次三次式 f と λ≠0 に対し

    P_f = Jac(f)(A,B,C) + λ W∂W ∧ (A∂A+B∂B+C∂C)

を使う。Jac(f) の規約は

    {A,B}=f_C,  {A,C}=−f_B,  {B,C}=f_A.

Jac(f) は Poisson であり、三変数 Euler 場はこの二次斉次 bivector を保つ。
従って P_f は全 f,λ について strict Jacobi を満たす。実 affine Pfaffian は

    Pf(P_f) = −3λ W f.

λ=-1,f=ABC は退化式 `3H W A B C` を持つ一般階数四の SNC toric 点である。
補助源 `G_m × Sym³(C³)` は既約なので、その像閉包を含む全 Poisson 成分を一つ選べる。
その成分は上の SNC 点を含む。採用済みの MPS SNC 成分定理により、この成分は
既知40成分の一つである。従って補助源の**全像**は既知40成分の有限閉和に入る。
これは `TWO_HYPERPLANES_CUBIC_SUPPORT.md` §2 の scalar-B 源にも用いた接続である。
個々の f の滑らかさを仮定せず、どの行政 ID になるかの追加同定も必要としない。

## 2. 二超平面を二節点接線へ衝突させる実族

g(A,C) を任意の二変数三次式、G(A,B,C) を任意の三変数三次式とし、

    f0 = −λABC + g(A,C),       fε = f0 + εG

と置く。ε≠0 に対して射影線形座標変更

    H=h,  A=x,  B=y,  C=h+εz,  W=x+εw

を行い、bivector を ε 倍する。座標行列は上記順序で行列式 ε² を持つ。
逆変換の z=(C−H)/ε,w=(W−A)/ε に現れる分母は最終的に全て消える。
実際、`a` を fε の A,B,C 微分の代入後の値とする規約で、全六括弧は

    {x,y} = ε f_C,
    {x,z} = −f_B,
    {x,w} = −λ(x+εw)x,
    {y,z} = f_A,
    {y,w} = −λ(x+εw)y + f_C,
    {z,w} = −λh w −λεz w − G_B(x,y,h+εz).

最後の式は

    (−λWC−f_B)/ε
      = [−λ(W−A)C−εG_B]/ε

という恒等式による。したがって全係数は ε と全源パラメータの多項式であり、
(h,x,y,z,w) の斉次二次式である。上記六括弧は任意の ε で strict Jacobi を満たす。
実 Pfaffian と射影退化式はそれぞれ

    −3λ(x+εw) fε(x,y,h+εz),
    −3λ h(x+εw) fε(x,y,h+εz).

ε=0 の bracket は、一般の f0 では h,x が節点の二本の接線となる配置を持つ。
重要なのは、六括弧自体が ε=0 でも定義された実族である点である。
分母を捨てた rational 同型だけ、接線だけ、又は有限 jet だけの議論ではない。

## 3. 既存の全六磁気源との係数恒等式

既存の `node_two_node_tangents` では h=Y−X,x=Y+X,y=Z なので、

    C_node(h,x,y) = hxy − (x−h)³/8.

従って λ=-1,g(A,C)=−(A−C)³/8 を選ぶ。
磁気二次式を

    R = r0 h²+r1 hx+r2 hy+r3 x²+r4 xy+r5 y²

とし、その三次原始式を

    G = 3r0 C²B + 3r1 CAB + (3/2)r2 CB²
        + 3r3 A²B + (3/2)r4 AB² + r5 B³

と置けば、厳密に G_B(x,y,h)=3R となる。
極限の fiber 座標を

    Z=−8z/3,    V=(z+w)/8

へ変更し、tensor を24倍する。この一定 GL2 変更の行列式は −1/3 であり、
いかなる r_i も分母に現れない。得られる全括弧は

    {x,y}=0,
    {x,Z}=64hx,
    {x,V}=−3hx+3x²,
    {y,Z}=24h²−48hx−64hy+24x²,
    {y,V}=3hy+6xy,
    {Z,V}=−64hV−3hZ+24R.

これは既存 report の `full_six_magnetic_source` と**全係数で一致**する。
実 affine Pfaffian は `−576x C_node` であり、全 r_i において非零である。
従って全対象点は一般階数四のまま保持される。

## 4. 閉包による全支持への結論と証明書

ε≠0 の全点は第1節の既知閉和を PGL と非零 scalar で移した点である。
各 Poisson 成分閉包は連結な PGL5 作用で保たれ、その有限和は閉である。
従ってこの**多項式実族**の ε=0 の全点も同じ閉和に入る。
全六 r_i とその全境界を残した係数恒等式により、対象源全体が被覆される。
既知成分への一般点所属から未知 fiber を推測する推論は用いていない。

exact producer は `code/classification_plane_cubic_toric_collision.py`、
有限契約は `jobs/plane_cubic_toric_collision.json`、
全式・入力 hash・producer bytes の証明書は
`reports/CLASSIFICATION_PLANE_CUBIC_TORIC_COLLISION_20261107.json`。
QQ 上で全 λ,g0,…,g3,r0,…,r5,ε を保ち、元模型と実族の全 strict Jacobi、
rational pushforward と多項式族の係数一致、全次数、Pfaffian、SNC toric 点、および
既存六源との全係数移送を確認した。実測3.221秒。

独立再構成はPLANE_CUBIC_TORIC_COLLISION_REVIEW.mdと同名REVIEW reportに保存した。
producerをimportせず、任意全十係数三次式の元模型、実逆座標行列、全十Jacobi、
全六磁気列の移送とPfを再構成し、SNC成分定理の適用も独立読解確認した。
残る三配置は `node_flex_tangent_through_flex`、
`cusp_flex_secant_through_flex`、`cusp_cusp_tangent_cusp_secant` である。
