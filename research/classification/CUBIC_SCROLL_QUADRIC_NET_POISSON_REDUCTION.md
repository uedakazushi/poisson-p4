# cubic scroll の quadric-net は Poisson 射：零 base の全排除

2026-11-07。主義務 G。純幾何 producer proof、別作者 review 前。
前提と全 P2-bundle は `CUBIC_SCROLL_QUADRIC_NET_GEOMETRY.md` に置く。
新しい全 G 次元上界・成分・E 帰属は主張しない。

## 1. 主張

smooth cubic scroll S と実既約被約 quintic D、一般 rank 二 restriction の全 source を考える。
`Y=Bl_S P4`、quadric-net morphism `p:Y→P2`、exceptional divisor E、strict transform D' を取る。

**必要還元。** p は P2 上のある projective Poisson tensor `tau` への Poisson morphism。
`tau` は恒等零でなく、その実 cubic 零因子は quadric-net discriminant line `L={c=0}` を含む。
従って base cubic は `c F2(a,b,c)` の形であり、F2 の可約・非被約・零係数を全て許す。
以下では全 symbolic-square と Euler gauge の矢印を証明する。

## 2. 全 raw lift を一つの Euler gauge で scroll-cone tangent にする

quadratic homogeneous lift Pi の projective tangency は、S 上で

`Pi#dq_i=lambda_i Eul`

を意味する。lambda_i は degree 二 homogeneous functions であり、Hilbert–Burch syzygies を満たす。
従ってこれらは normal bundle N の一つの global section lambda に等しい。
Euler gauge `Pi→Pi+Eul∧W`、W は ambient linear vector field、とすると
lambda は W の normal projection だけ変化する。

`H1(S,T_S)=0` を直接確認する。
`rho:F1→P1` に対して relative tangent sequence は

`0→O(2C0+f)→T_S→rho*O(2)→0`。

`rho_*O(2C0+f)=O(1)⊕O⊕O(−1)`、R1rho_* は零なので左の H1 は零。
右の H1 も零であり、従って H1(T_S)=0。
normal sequence から `H0(T_P4|S)→H0(N)` は全射。
restricted Euler sequence と H1(O_S)=0、H0(O_S(H))=C5 により
`H0(T_P4)→H0(T_P4|S)` も全射（ともに dimension24）。

従って一つの ambient W で全 lambda_i を同時に消せる。
この gauge において raw polynomial vector `Pi#dq_i` の全成分は I_S に入る。
projective degree 零 functions の bracket は Euler gauge で変わらないので、後の net map は同じである。
この lift に canonical divergence-zero 条件を同時に課してはいけない。

## 3. 全 degree 四 symbolic square は q の二次式

intrinsic conormal bracket の可換性は、chart x0 で

`x0{q_i,q_j}−2q_j{q_i,x0}+2q_i{q_j,x0}∈I_S^(2)`。

§2 の gauge では後二項は I_S² に入るので、
`{q_i,q_j}∈H0(I_S^(2)(4))`。

blowup の全 P2-bundle 表示では `4H−2E=2B`。
S は smooth なので、degree 四 polynomial が S の全 generic normal direction に order 二以上で消えることは
blowup 上の `4H−2E` の section に正確に対応する。
従って

`H0(I_S^(2)(4))=H0(Y,O(2B))=H0(P2,O(2))`。

p の全 fibers が connected P2 なので p_*O_Y=O_P2。
右辺の六 basis は正確に q_iq_j。よって全 coefficients を保持して

`{q_i,q_j}=F_ij(q0,q1,q2)`、F_ij は constant-coefficient quadratic。       (1)

chosen parameter chart の範囲だけの同一視ではない。
S の codimension 二の bad locus や square の saturation を無断で捨てていない。

## 4. 全 morphism 上の Poisson 性

(1) は P2 の quadratic homogeneous bivector lift を定める。
P2 上では任意の bivector は自動的に Poisson であり、これを tau とする。
q_i/q_j の全 degree 零 ratios の bracket は (1) の quotient rule により tau の bracket と一致する。
従って P4 の rational net map は Poisson。
S に沿う intrinsic conormal 可換性により Pi は blowup Y へ regular に持ち上がり、
二つの regular tensor が dense open で一致するので p は Y 全体の Poisson morphism。
E と D' はそれぞれ Poisson divisors である。
E は blowup での直接構成、D' は generic symplectic Poisson の反標準因子の irreducible component から従う。

## 5. tau=0 は全 source で不可能

tau=0 と仮定する。base の任意の local coordinates f,g の Hamiltonian vector fields
`X_f,X_g` は vertical であり、Poisson Jacobi により可換。
generic fiber F=P2 上に regular vector fields として制限できる。
両者は Poisson divisors E,D' を保つので、F 上の smooth conic `E_F` と line `D'_F` を保つ。
generic base point を L と hyperplane section の全-zero locus の外に取れば、この conic と line は確かに存在する。

smooth conic の preserving group は PGL2。
line の conic への intersection が二相異点なら pair の Lie algebra は一次元 torus。
double point なら二次元 Borel algebra で、その非零元の centralizer は常に一次元。
従ってどちらの全場合でも、pair を保つ commuting vector fields 二つは一次元以下の span しか持たない。

一方 generic fiber の E_F∪D'_F の外では Y の Poisson rank は四。
df,dg は独立なので X_f,X_g も独立でなければならない。矛盾。
従って tau≠0。
line が conic に tangent である全境界もこの議論に含み、secant open だけを調べていない。

## 6. conic-bundle discriminant L は base cubic の一成分

p|E は Poisson morphism。
各 point y∈L の singular conic fiber は二本の異なる lines の和で、node point が一つある。
total E は smooth で、その node point では differential `d(p|E)` の rank は一。
局所 model `uv=t` で dt=u dv+v du となることからも直接分かる。
もし tau(y) が rank 二なら、この rank 一 differential で pi|E の bivector を tau(y) へ送れず矛盾する。
従って L 上で tau は恒等零。
tau≠0 の anticanonical cubic は `c F2` と因子分解する。

## 7. 残る数学

base quadratic F2 は irreducible conic、split/repeated lines、c 自身の重複を全て許す。
relative hyperplane D' には長さ四の必須 full-P2 fibers と追加 isolated zero の補集合がある。
そのため本還元から flat P1-bundle、SNC 全域、全 residue normal form 又は dimension≤24 を推測しない。
残る A=(2,5) 二型を解くには、この全 Poisson P2-bundle の connection と全 full-fiber boundary を扱う必要がある。
cubic cone はこの morphism/normal-bundle argument と別の義務である。
