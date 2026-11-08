# nil31短鎖 r=2/3：二parameter modular plane 全体の flag 障害

2026-11-08、PRIMARY E、`/root/quartic_review` の純数学命題。
未決全源 `MODULAR_NIL31_SHORT_SUPPORT.md` §7 の所属は**未解決のまま**である。
新 CAS、source/operator、profile、jet、minor 又は数値点は実行しない。
以下は一つの固定 modular 方向だけでなく、明示した全二parameter plane の
任意の ramified regular curve と全 quadratic Poisson 源を排除する。
任意の modular 変形、任意 MPS 族、既知成分所属を否定する命題ではない。

## 1. 全対象と結論

canonical source の全 t,d について

    Pxy=d z−t y+2t zw−y²/3, Pxz=t z−2yz/3,
    Pxw=d+yw/3, Pyz=−2z²/3,
    Pyw=(y+zw)/3, Pzw=2z/3                         (1)

と書く。div P=N0=y∂x+z∂y+∂w、rank四条件は t!=0。
元六parameter源から (1) への determinant 一 affine inverse は既証明のまま保持する。
td!=0 は一つの既約22次元軌道、d0 はその明示多項式源の閉集合である。
その E被覆には全源を捨てない actual curve 又は別の認識が依然必要。

R=C[[s]]（有限 ramified extension も許す）、K=Frac(R) とする。
rho,sigma∈sR、sigma!=0 とし、次の**全二parameter** affine modular field を取る。

    N_(rho,sigma)=y∂x+(z+rho*w+sigma*x)∂y+∂w。       (2)

命題：次数二以下の strict affine Poisson bivector P_s が R 上 regular、
P_0=(1)、div P_s=(2) となることは無い。
特に旧候補 rho=epsilon,sigma=epsilon² の全 quadratic 源を排除する。
旧 `NIL31_SHORT_SEMISIMPLE_CANDIDATE.md` の fixed three-variable pi の不整合を
任意 pi へ外挿したのではなく、以下の全源 argument を新たに用いる。

sigma零の全 curve、(2) 以外の matrix/translation deformations、
異なる h-Casimir chart を要する全親族は本命題の補集合である。
これを未確認の Slodowy normalform の完全分類とは呼ばない。

## 2. 全 BV 不変 quadratic 源を一度に扱う

Poisson から BV identity δ[P,P]=2[δP,P] により L_N P=0。
sigma の平方根 lambda を有限拡大で取り、B=z+rho*w+sigma*x として

    a=B+lambda*y+rho/lambda,
    b=B−lambda*y−rho/lambda,
    z=z, w=w

を取る。determinant は −2lambda³、K 上で可逆。
inverse は y=(a−b−2rho/lambda)/(2lambda)、
x=(a+b−2z−2rho*w)/(2sigma)。
これは各 generic fiber の actual affine change であり、s0で可逆とは仮定しない。
divergence は constant Jacobian の affine change で自然に移送され、

    N=lambda*a∂a−lambda*b∂b+∂w。

quadratic polynomial tensors を diagonal weight に分ける。
nonzero weight では constant nonzero scalar+∂w が有限 nilpotent geometric series で可逆。
weight零では∂w不変でなければならない。
従って全 L_N P=0 polynomial 源は w 非依存、weight零である。
特に以下の三つの全括弧は正確に

    Paw=a D(z), Pbw=b E(z), Pzw=f(z)+k ab、             (3)

D/E は全 affine、f は全 quadratic、k は任意 scalar in K。
他の全括弧を零に置いていない。
例えば全 source はさらに Pab=g(z)+C ab、Paz=a A(z)、Pbz=b B(z) を持つ。
div 条件は B=−A、C=−lambda−A' を与えるが、**A' は自由**である。
A' を零と扱って初期 quadratic Pyz を失う誤った shortcut は使わない。
元 Jacobi の追加必要式も本障害には不要。

div の w 成分は generic coordinates で厳密に

    D(z)+E(z)+f'(z)=1。                              (4)

全 constant/linear/quadratic coefficients、全 affine particulars、全 pole orders を残した恒等式である。

## 3. 元係数に戻した分母無しの必要恒等式

元原点 x=y=z=w0 では a=rho/lambda,b=−rho/lambda。
inverse を使い

    constant(P_old^yw)=rho*(D(0)+E(0))/(2sigma)。

また ab=B²−(lambda*y+rho/lambda)² は元 z の一次単項式を持たない。
従って coefficient_z(P_old^zw)=f'(0)。
rho!=0 のとき、全源に

    (2sigma/rho)*constant(P_old^yw)
      +coefficient_z(P_old^zw)=1                    (5)

が成り立つ。分母を払えば

    2sigma*constant(P_old^yw)
      +rho*coefficient_z(P_old^zw)=rho。

定数項は元 affine 原点での値、coefficient_z は元 quadratic bracket の z一次係数。
保存 coordinate pole を勝手に打ち切った endpoint matching ではない。
旧 rho=epsilon,sigma=epsilon² なら (5) は
2epsilon*constant(Pyw)+coefficient_z(Pzw)=1。
regular endpoint では左の第一項が消え、(1) の第二項2/3と矛盾する。

## 4. 全相対 valuation と zero complement

BV 不変 source は generic coordinates で w 非依存なので、元に戻すと

    L_(c0)P=0, c0=∂w−(rho/sigma)∂x。                 (6)

この translation は linear semisimple part と commute する zero-eigen component。
全 tensor が polynomial なので、N 不変から個別 translation 不変へ移す §2 の weight argument が有効。
未確認の Jordan 分解不変性を仮定していない。

target (1) の constant translation stabilizer は正確に C∂x。
実際 V=alpha∂x+beta∂y+gamma∂z+delta∂w に対し、
L_V Pyz=−4gamma*z/3 から gamma0、
L_V Pxz=−2beta*z/3 から beta0、
L_V Pyw=delta*z/3 から delta0。alpha は任意。

rho0 なら (6) は∂w不変で endpoint もそうなるため即矛盾。
rho!=0 で rho/sigma∈R なら (6) を s0へ特殊化し、
∂w−(rho/sigma)(0)∂x が target の symmetry となるので矛盾。
残る whole complement は rho/sigma が pole を持つ場合のみ。
このとき sigma/rho∈sR だから、regular P の (5) を s0へ特殊化すると

    coefficient_z(P0^zw)=1。

しかし (1) は2/3。矛盾する。
全 rho/sigma の相対 valuation、任意有限 ramification、全 quadratic coefficients を一括で尽くした。
scalar unit kappa(s) 倍の modular field を許しても右辺が kappa に変わるだけで、
nonzero scaled target の2/3比は一致しない。

## 5. 正確な進展と次の義務

この全 modular plane の中では、translation flag を target の∂xへ潰す必要がある一方、
その flag degeneration が (5) の residue coefficient を1へ強制する。
target の2/3との不一致が、fixed candidate を更に高次へ延長しても直らない構造的原因である。
正しい actual curve は modular matrix/translation の形をこの plane の外へ動かすか、
別の既知全 source を直接認識する必要がある。
同じ plane の新 profile、別の相対次数、fixed three-variable tensor の再試行はしない。
全 shortnil r2/3 E所属、新成分性、任意 MPS 変形の存在/不存在は未解決として保持する。
