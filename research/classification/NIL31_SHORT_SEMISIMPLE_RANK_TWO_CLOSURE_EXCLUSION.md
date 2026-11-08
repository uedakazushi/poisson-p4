# nil31短鎖 r=2/3 と marked 重複超平面・rank二半単純 modular 層の閉包

2026-11-08、PRIMARY E、`/root/quartic_review` の純証明。
新 CAS、数値実行、source/profile、operator、minor、jet 又は GB は無い。
本稿は既知閉和への正の帰属ではなく、失敗した候補の正確な幾何的原因を全源へ拡張する。
全 r2/3 E所属は未解決として保持する。

## 1. 正確な命題と入力

target は `MODULAR_NIL31_SHORT_SUPPORT.md` §7 の全 canonical source(15)、t!=0、d任意。
固定 affine chart h1 で target の div は N0=y∂x+z∂y+∂w。
先行純証明 `NIL31_SHORT_RANK_TWO_MODULAR_FLAG_OBSTRUCTION.md` は

    N_(rho,sigma)=y∂x+(z+rho*w+sigma*x)∂y+∂w,
    rho,sigma→0, sigma generic nonzero

の**全二parameter plane**に対し、全 regular quadratic Poisson endpoint を排除した。
使った証明は全 BV invariant source の係数恒等式と translation flag の二分岐であり、
fixed three-variable pi の旧負証明を全 pi へ外挿したものではない。

今回の命題：target への任意の regular quadratic Poisson curve について、
generic modular affine linear part A が semisimple rank二となることは無い。
任意 c を保持し、c∈im A と c∉im A を分ける。
従って actual multiple hyperplane を標識したこの rank二半単純全層の閉包にも target は入らない。
chosen infinity hyperplane が単純な一般 SS2 source の unmarked global closure へは拡張しない。
他の rank三/四 modular 層、非半単純 mixed 層、異なる実因子型の親全族は本排除の範囲外。

## 2. 任意の rank二半単純 affine 場の whole regular normalform

R=C[[s]]、K=Frac(R)、有限 ramification を許す。
N=c(s)+A(s)x が R 上 regular、N0 が上記 target、generic A は semisimple rank二とする。
div(div P)=0 から tr A=0。従って generic eigenvalues は lambda,−lambda,0,0。
sigma=tr(A²)/2∈sR は generic 非零で、A³=sigma A。

N を affine linear functions V*=R·1+R·x+R·y+R·z+R·w 上の R-linear endomorphism と読む。
generic c∈im A なら、translation で c を零へ移せるから rank N=2。
target は rank N0=3 なので rank≤2 の閉条件に反し、この全補集合は不可能。
以下 generic c∉im A とする。
generic N は二非零 eigenlines と zero J2+J1 の直和で、rank N=3。
target も rank3、従って DVR 上で N の rank は一定。
nonzero 3minor が unit、4minors が零なので image と kernel はともに direct-summand free bundles。
これは仮定した normalization-flatness ではなく、元五次元 linear map の constant rank の帰結。

generic では 1∈im N である。
実際 c∉im A なら ker A* 内に ell(c)!=0 の covector があり、尺度を変えると N ell=1。
image の locally-free quotient で section1 の class は generic零、従って R 上も零。
ゆえに N w_s=1 を満たす regular affine function w_s を取れる。
special fiber で w_s−旧w∈ker N0。
kernel bundle の reduction は ker N0へ全射なので、その差を kernel section で除き
**w_s→旧w**としてよい。

さらに generic affine functions 上で

    N⁴=sigma N²,
    im(N³−sigma N)⊂R·1。                              (1)

後者は非零 eigenlines では零、zero J2 の上では −sigma N、J1上では零なので従う。
両者とも元 R-linear matrices の恒等式だから全 Rへ延長する。
旧 x を固定し

    y_s=Nx,
    rho=(N³−sigma N)x∈R·1,
    z_s=N²x−sigma*x−rho*w_s

と定義する。rho0=0、N z_s=0、N y_s=z_s+sigma*x+rho*w_s。
(x,y_s,z_s,w_s) は s0で正確に旧(x,y,z,w)なので determinant が unit。
従ってこれは**regular actual affine coordinate change with identity endpoint**であり、
全 P の quadratic degree、regularity、strict Jacobi と divergence を保つ。
この座標で N は厳密に先行全 plane の N_(rho,sigma) に入る。
rho/sigma の全相対 valuation も先行命題で保持されている。

この証明は一つの Slodowy sample、fixed matrix candidate 又は generic basis の極限を使わない。
任意の rank二半単純+c regular deformation を全源のままこの plane へ戻す。

## 3. target 固有の flag/residue 矛盾

normalform へ戻した P の endpoint は変わらない。
先行 fullsource coefficient identity は rho!=0 で

    (2sigma/rho)*constant(Pyw)+coefficient_z(Pzw)=1。    (2)

generic translation symmetry は c0=∂w−(rho/sigma)∂x。
target の constant translation stabilizer は正確に C∂x。
従って rho0 又は rho/sigma regular は c0 の非零∂w成分により矛盾。
残る全場合は rho/sigma が pole を持ち、sigma/rho→0。
(2) は regular endpoint に coefficient_z(P0zw)=1 を強制するが target は2/3。
これで任意の generic rank二半単純 affine modular curve が排除された。

既済 r1/3 主枝の正しい semisimple deformation とは矛盾しない。
その別 canonical tensor は w 非依存であり、target の translation flag が異なる。
同じ nilpotent modular matrix であっても、今回の二つの tensor-specific 条件
「constant symmetry は∂xのみ」と「Pzw の z係数2/3」を抜いて適用してはならない。

## 4. marked multiple hyperplane の global closure への橋

ここで global closure を述べる層は、actual multiple hyperplane H を標識したものだけである。
既存 `MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` の冒頭と§2は、
明示的に P=C+L+Q、**Pf(Q)=0** を仮定した全 rank二半単純源を被覆する。
h-Casimir quadratic representative の actualD は
h·Pf(h²C+hL+Q) で、h0 の係数 Pf(Q) が零だから実際に h² が割る。
したがってこの既存 theorem の chosenH は確かに multiplicity≥2である。
これは semisimple rank二から自動的に従うと仮定したものではない。
projective point の curve を仮定した場合、hyperplane と D の incidence H²|D は closed、
dual projective space は properなので finite base change 後に H_s を延長できる。
target の actualD は非零定数倍の h²z(hy−zw)、唯一の repeated hyperplane は h。
従って H0=h、regular PGL change（endpoint identity）で H_s を h に固定できる。

この fixed h chart の quadratic h-Casimir representative も regularに取れる。
標準 divergence-free homogeneous lift Pi は係数の linear operation で regular。
Pi#dD=0 と factor divisibility から Pi#dh=hX、h は coordinate unitなので X も R 上 regular。
Euler correction Pi−E∧X の h-column は零、全 remaining quadratic coefficients は regular。
h1 に制限すると次数二 strict affine Poisson familyである。
これが前節の全源 normalform/contradiction の対象になる。
scalar frame を endpoint 非零係数で正規化すれば target は(15)そのもの。

constructible algebraic stratum の closure へ入るなら、標準的な curve/DVR criterion と
上の proper marked-factor incidence によりこのような generic-in-stratum curve を選べる。
従って **shortnil r2/3 source は marked 重複H・rank二半単純 modular 全層の閉包外**。
canonical td!=0 一軌道だけでなく、t!=0,d0 を含む target の二条件により同じ結論となる。
generic source の係数を捨てず、hyperplane の変動も unit coordinate change で明示的に戻した。

## 5. E課題に対する正確な制限

「generic modular が既済 rank二半単純+c 型へ入る actual curve を作る」という
従来の正の提案は、この target に対して全源で不可能だった。
fixed pi 又は単一 flag だけの計算失敗より強い構造的障害である。
この層へ向かう別次数・別 matrix normalform・higher jets を追加しない。
実帰属を得る次手は、rank四半単純/別 mixed 型、repeated hyperplane が消える親全族、
又は既知成分全 source の異なる actual recognition でなければならない。
generic chosenH が単純な SS2源では、標識Hが target の z に近付く可能性を排除していない。
任意 MPS 閉和からの非所属、全 r2/3 source の新成分性、全 E完成は未主張。
