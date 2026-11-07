# 三つの high-m 残余：κ非零の actual normalization と全 Jacobi の三角還元

2026-11-07。主義務 G、producer `plane_collision_review`。
対象は `PLANE_RANK_ZERO_RANK_TWO_HIGH_MULTIPLICITY_BOUND.md` の最後の三つ、
Borel \(M=bH+cE\)、\(U=\partial_a\)、\(T=\rho b\partial_a\)、
\(\rho=1,5/3,3\) の **全 higher source** である。
本ノートは旧 flag rowspace criterion を multiplier degree 上げで繰り返すものではない。
実際の PGL transformation で symplectic anchor を固定し、その全源の affine Jacobi を扱う新しい幾何還元である。

## 1. 既受理の全係数と、本質的 open の完全性

三つの P0 は

\[
 \{a,z\}_0=(\rho+1)bz+cw,\qquad\{a,w\}_0=(\rho-1)bw
\]

だけが非零。全 P1 kernel は14,14,15係数、P2 は全九係数を持つ。
d を P1 の \(\{a,w\}\) における \(z^2\) 係数、r を P2 の \(\{b,c\}\) における \(z^2\) 係数とする。
保存した original350-row functional は、任意の全 Jacobi solution 上で

\[
 D|_{w=0}=\kappa z^5,\qquad\kappa=-dr
\]

を与える。後式は κ の rowspace remainder が正確に source pair (8,20)/(8,21) の負号であることによる。
κ0 なら w が実 quintic D の因子となり、既約被約に反する。
従って今回対象の全実 rank 四・既約 D 枝は \(d\ne0,r\ne0\) に入る。
κ0の補集合を未処理として落としてはいない。

chart z1,w0 では Hamiltonian \(X_w=d\partial_a\)、
ρ1,5/3の他の modw 係数の零性は追加 P1係数ℓについて dℓ0を与え、d非零からℓ0。
従って全三型の coisotropic quotient の (b,c) bracket は r となる。
これは有限 sample の Darboux point を選ぶものではなく、全 source の二つの非零係数である。

## 2. 全 d,r に有効な実可逆 PGL transformation

次の parameter-dependent diagonal coordinate change と tensor scalar を用いる。

\[
 a'=\frac r{d^3}a,\qquad b'=\frac1d b,\qquad c'=\frac1d c,
 \quad z'=z,\quad w'=w,\qquad\sigma=\frac{d^2}r. \tag{1}
\]

全係数は \(\mathbb C[d^{\pm1},r^{\pm1}]\) 上で正則、行列は可逆、逆は対角因子と σ を逆にするだけ。
P0 は literal に保存される。
一般に \(a'=Aa,b'=Bb,c'=Bc\)、\(\sigma=B/A\) は P0 を保ち、
\(d'=Bd\)、\(r'=B^3r/A\) と作用する。
\(B=1/d,A=r/d^3\) が正確に \(d'=r'=1\) を与える。
grading を保つので全30/9 source を保ち、P1 の complete kernel も Jacobi の共変性から保つ。
他の unknown coefficient を零にしたり、Poisson family を一点へ特殊化したりしていない。

全係数の検査は二整数の weight だけでも行える。
座標の (d,r) weight は
\(a=(-3,1),b=c=(-1,0),z=w=(0,0)\)、scalar は (2,−1)。
raw bracket \(\Pi^{ij}\) の monomial \(x^e\) の係数は

\[
 (2,-1)+s_i+s_j-\sum_k e_ks_k
\]

の weight で変換される。
P0 の全非零係数は weight零。
common-U の十二 P1係数は weight (−1,0)、追加 P1係数は (1,−1)。
P2 の \(\partial_a\wedge\partial_b/\partial_c\) 六係数は (−2,0)、
\(\partial_b\wedge\partial_c\) 三係数は (0,−1)。
全150係数でこの homogeneous weight 一致を検査出来るので、
各 normalized coefficient と inverse coefficient の全 map が Laurent monomial で与えられる。

従って実 Poisson source の d,r非零 open は、
\((\mathbb C^*)^2\) と **d=r1における全 normalized-source Jacobi locus** の積と同型である。
二つの consumed parameter を後の dimension 上界へ必ず戻す。

## 3. 低次数を保った affine Gaussian source 還元

保存済みの full350 Jacobi coefficient polynomials へ d=r1を literal に代入する。
残る unknown は21,21,22で、全 polynomial は次数高々二。
全 equation の QQ-linear rowspace 内にある affine-linear equations を取り、
nonzero rational pivot で未知係数を解いて全 equation へ代入する。
各 solved equation について、その時点の **全 original350 equations の定数 row functional** を保持する。
これは source を動かす actual affine coordinate elimination で、zero locus と residual source の間の全 map/inverse は明示的である。

必要ならこの affine elimination を繰り返すが、毎回少なくとも一 unknown が消え、
回数の先験上限は22、polynomial degree は常に二以下。
ideal multiplier の次数、normal jet の次数、GB の範囲を増やす操作は無い。
各段階で全350 residual equations を残し、linear consequence のみで全源を含む有限 affine space へ埋め込む。

normalized source の ambient dimension が十以下になれば、それ以上の equation を解かずに止めてよい。
全 normalized Poisson locus はその residual zero locus のため dimension も十以下。
affine equation が無く unknown 十超なら exact residual を保存して停止する。
constant 非零の consequence が出れば source は空である。
この bounded 手続は full nonlinear Jacobi の全分類や非被約構造を求めるものではない。

## 4. 全 G 上界へ戻す正確な数え方

各 fixedρ の全 first-symbol affine source dimension は七。
plane choice六、tensor projectivization一を戻すと offset は十二。
§2 の二 anchor parameter は normalization で消費したため、
normalized residual source dimension 上界を t とすれば全像は

\[
 7+(t+2)+6-1=14+t
\]

以下である。従って t≤10 が最弱の十分条件である。
一点の ambient tangent rank や generic orbit から全 family を推論していない。
κ0は既約 D から排除済み、κ非零全体は actual (1) と逆によって復元される。
新しい有限証明書の成功前には、この三源の G完了を主張しない。
