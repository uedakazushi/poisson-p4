# Ga9 linear-w：既知の親候補と二次方向の幾何

2026-11-07、PRIMARY E。対象は `GA9_LINEAR_W_QUADRATIC_FRONTIER.md` の全23軌道源である。
全六normalの二次支持は二平面 P と直線 L の和で、射影局所次元の上界25だけが既済である。
本稿は既存証明書の適用範囲と手で展開できる不変量を整理し、
高次jetを再開せず次の actual whole-family 認識の候補を絞る。
新たな Poisson 曲線や成分性は主張しない。

## T39/T40 の既存 negative はどの chart 上の判定か

入力 `CLASSIFICATION_HYPERPLANE_QUARTIC_NIL211_ZERO_C_LINEAR_W_GRAPH_20261107.json`
の SHA256 は `e97e517050777ef0b5e318fcb73146210d55fce84d78395fc331f68abdc3d0da`。
この全源の密軌道の一点を、affine 点 (0,1,0,0) で評価している。
保存された定数行列は

    K=[[0,3/4,8,4],[-3/4,0,1/4,0],[-8,-1/4,0,1],[-4,0,-1,0]]

で、Pf(K)=7/4、det(K)=49/16 である。
従ってこの評価は非退化な開集合内にある。
translation の determinant は1、保存 Darboux frame の determinant は4/7であり、
perspective parameter は (0,−1/6,0,0) である。
perspective は固定可逆 projective 行列で、βに依存する逆元を導入していない。
全61座標への再構成と元134式は旧reportに保存済みである。

`Q51_UNIVERSAL_SECTION_20261101.md` と `Q52_UNIVERSAL_SECTION_20261101.md` の既済定理は、
T39、T40 の**被約 S61 閉包そのもの**を A16 上の多項式 graph と同定している。
有理係数の定数分母はあるが、βのある多項式を反転した graph chart ではない。
この既済の全閉包同定と S61 normalization の同値性を入力とすれば、
旧検査の第3座標（index2）の残差

    T39: −3/112,   T40: 99/112

は chart 外を見た negative ではない。
従って、この全域 graph 同定の適用範囲において canonical 軌道は両成分に含まれない。
本稿で全134式や graph を新しく計算したわけではなく、既存の同値性定理を再証明したわけでもない。

## 五つの24次元 MPS 成分の閉因子分割

全sourceの homogeneous bivector を B、Euler vector を E とする。
退化五次は `D=coefficient(E∧B∧B/2)` であり、全五座標で

    D=Σ_i (−1)^i x_i Pf(B の第i行列を除いた4×4行列)

と計算する。同じ scalar convention を全てに使う。
以下は公開 MPS source 式からの直接展開で、新しい一般点の標本検査ではない。

T02、T03、T16、T21 の式は `MPS_S61_RAW_CATALOGUE_20261103_v62` の同ID JSON、
T17 は `MPS_SMOOTHING_COMPLETION_20261103_v62/T17.json` にある。
その period を p_i、smoothing を a=t0、b=t1、c=t2 と書く。

T02では λ=p0+(p2+p3)/2 とおくと

    D=λ x2 x3 x4 (p1 x0 x1+a x4²).

T03では λ=2p0+p2+p3 とおくと

    D=λ x2 x3 x4 ((p1−p2)x0 x1/2+a x3 x4).

T16では μ=(p2−p1)/2、λ=p0+(p1+p2)/2 とおくと

    D=x4 (μ x0 x1+a x4²)(λ x2 x3+b x4²).

T17の completed source には p1 の逆元を含む correction があり、全体で

    D=x4 (p1 x0 x1+a x4²)
         ((p0+p2)x2 x3+b x0 x1+ab x4²/p1).

T21では

    D=p0 x3 x4 (p1 x0 x1 x2+c x0³−b x1³+a x2³).

従って全一般源の因子分割は順に3H+Q、3H+Q、H+Q+Q、H+Q+Q、2H+Cである。
分母を持つ T17 は p1≠0 の既済源を使うが、この結論は全成分閉包にも延びる。
実際、所定次数の因子を持つ五次形式は射影空間の積からの積写像の像であり閉である。
H+既約quartic はこれらの閉像のいずれにも属さない。
ゆえに Ga9 の密な23軌道はこの五つの24次元成分に入らない。

残る既知24次元 MPS 親候補は T36、T37、T38 である。
これらの既存 plane normal-cone の部分族だけで非所属とは結論しない。
103番目の24次元成分との区別は `GA9_LINEAR_W_MARKED_PAIR_OBSTRUCTION.md` の
marked-pair 命題で別に処理する。

## 二次の二平面は H 保持、孤立直線は H を破る一次方向

ここでは新しい Poisson 曲線を作らず、保存された full150 一次方向と Pfaffian の微分だけを使う。
二次証明書の整数scale4の点は h=0 上で

    Pxy=−y², Pxz=R+2yz, Pxw=−4wy,
    Pyz=Pyw=0, Pzw=4y²,
    R=44w²+28wy+20y².

退化五次の h=0 における一次変動を、h 行の変動 α_i=δP_hi で書くと

    δD|h=0=4y³ α_x+3wy² α_z+α_w*y*R+3α_w*y²z.

旧 normal e0 の h 行は

    α_x=16hy−3yw,  α_z=−8h²−3hw+3y²,

他は零で、δD|h=0=−3y⁴w となる。
canonical F|h=0=−4y⁴ なので、これは H の一次変動を (3/4)w と取って吸収できる。
normal e1 の h 行は α_x=−272h²−12hw、α_z=4hy、他は零で、δD|h=0=0。
従って全二平面 P=span(e0,e1) は H+quartic incidence の一次接空間内にある。

一方、孤立直線 L の保存した全150方向では、h=0 の行に

    α_x=−2702133*y²/5+3969*yz+18522*yw+72765*w²,
    α_y=0,
    α_z=37044*y²+16632*yw,
    α_w=−7938*y²

がある。δD|h=0 の y³w² 係数は

    4*72765+3*16632+44*(−7938)=−8316≠0.

従って δD|h=0 は y⁴ で割れず、H の一次変動で吸収できない。
この方向は simple-H incidence を一次で離れる。
これは quintic 親成分への候補方向を区別するだけで、実 Poisson 曲線への延長を主張しない。
真の局所成分が必ずこの方向を持つという主張でもない。

## 残す義務

既済の二次支持 P∪L は m-primary でなく、初期支持から新成分性を証明する経路は止めた。
T36/T37/T38 の whole closed source に対する actual 認識、又は
この点を通る未知の全既約 source の構造が残る。
fixed-N orbit-limit、commuting twist、D5係数変更、三次以降の jet を自動再開しない。
本稿の因子恒等式・一次 incidence 計算は手で検算した内容で、独立算術再現とは区別する。
