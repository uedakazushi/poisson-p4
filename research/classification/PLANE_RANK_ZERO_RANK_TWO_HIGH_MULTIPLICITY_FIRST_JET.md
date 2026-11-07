# 被約実 D・rank 零 plane・rank M 二：高重複度の全 first symbol

2026-11-07。主義務 G、producer `plane_collision_review`。
実退化因子 D が既約被約 quintic、\(S=(z=w=0)\) が rank 零 plane とする。
\(S\) 沿いの D の一般位数 \(m\ge3\)、co-Higgs 型 k=1、
係数写像 \(M:\mathbb C^3\to\mathfrak{sl}_2\) の rank が二の **全 first symbol** を扱う。
本ノートは幾何による全源還元であり、higher-weight fiber の上界は別の有限証明書に分ける。
rank M 三・一以下、k=0 は明示補集合である。

## 1. 全 m での可換 conormal と共通 Lax 式

新しい既受理命題 `REDUCED_DEGENERACY_ZERO_SURFACE_ABELIAN_CONORMAL.md` と別作者監査は、
被約実 D の rank 零曲面で \(\{I,I\}\subset I^2\) を全重複度について与える。
従って旧 m二入口の source 35/30/9 と k=0/1 共通線 subsheaf 分岐を、今回もそのまま使用出来る。
旧 m二で除いた leading-Pfaffian 零の点・族は全て復活させる。

\(U\) を O(1) inclusion の非零 constant vector、\(p=[U]\)、
\(n=(z,w)\)、\(E_n=z\partial_z+w\partial_w\)、\(q_A(n)=\det(n,An)\) とする。
trace-free linear lift T と \(S\in\mathfrak{sl}_2\) に対し全 weight 零 symbol は

\[
 \Pi_0=T\wedge E_n+U\wedge(Mn\cdot\partial_n)+q_S\partial_z\wedge\partial_w.
\]

\(M_0=M(U)\)、\(R=S+M_0/2\) とすると、全 projective weight 零 Jacobi は正確に

\[
 T(M)-tM+[R,M]=0,\quad[T,U]=-tU,\quad[R,M_0]=0. \tag{1}
\]

これは旧 `PLANE_RANK_ZERO_MULTIPLICITY_TWO_FIRST_JET.md` §2 で導かれた rank・m 非依存の式である。
normal order 二の actual Pfaffian は非零 universal scalar を除き

\[
 D_2=q_M(n)\det(E_{\rm base},T,U). \tag{2}
\]

実際、\(\Pi_0^2/2\) で残るのは \(T\wedge E_n\wedge U\wedge Mn\) だけ。
rank M 二なので \(q_M\not\equiv0\)、従って m≥3 は (2) の determinant の恒等零と同値である。
以下はこれを全 Lax source 上で課し、全 zero/row/resonance complement を保持する。

W=im M の二つの全 AdPGL2 型は
\(W_N=\langle F,E\rangle\)、\(W_B=\langle H,E\rangle\)。
\(H=\operatorname{diag}(1,-1)\)、\(E=\left(\begin{smallmatrix}0&1\\0&0\end{smallmatrix}\right)\)、
\(F=\left(\begin{smallmatrix}0&0\\1&0\end{smallmatrix}\right)\)。
Killing orthogonal line が nonisotropic/isotropic の全場合に対応し、それぞれの normalizer は CH、\(\langle H,E\rangle\)。

## 2. \(M(p)\ne0\)：全四 high-m source

旧 nonkernel first-jet note の正規化前の全 s,ℓ,τ を保持した式へ (2) を課す。
nondegenerate W の M0 semisimple/nilpotent の両場合では determinant の独立係数から
\(s=\ell=0\) が強制される。
Borel、M0 semisimple では \(\Delta=4\tau-3s=0\)、\(\ell=0\)。
Borel、M0 nilpotent では \(s=\ell=0\)。
例えば前二つはそれぞれ \(s a(b-c),\ell(b-c)^2\)、\(sac,\ell c^2\) の独立係数、
Borel semisimple は \(\Delta ac,\ell c^2\)、nilpotent は \(sab,\ell b^2\) から得る。
係数が zero になる任意の τ を捨てていない。

全 symbol は次の四つになる。K は raw \(q_S\) を表す。

| 全 source | U | M | T | K |
|---|---|---|---|---|
| NDSS0 | \(\partial_b+\partial_c\) | \(bF+cE\) | 0 | \(-(z^2-w^2)/2\) |
| NDN0 | \(\partial_b\) | \(bF+cE\) | 0 | \(-z^2/2\) |
| BSS0(τ) | \(\partial_b\) | \(bH+cE\) | \((-2\tau a,4\tau b,-2\tau c)/3\) | \((1-2\tau)zw\) |
| BN0(τ) | \(\partial_c\) | \(bH+cE\) | \(2\tau b\partial_c\) | \((1/2-\tau)w^2\) |

τ は全 C。旧 s≠0 や ℓ≠0 を使う shear/scalar normalization はここへ流用しない。
先行 whole image-plane/M0 normalization は base/normal の可逆変換と非零 tensor scalar だけであり、
全 zero/nonzero complement と逆変換は旧 note §§1–4 で保持済み。
上の high-m 条件を課した後に未知 higher coefficients を固定していない。

## 3. \(M(p)=0\)：high-m は正確に R=0、row は二係数全て

\(U=\partial_a\)、\(M=bM_1+cM_2\)、W=im M とする。
全 Lax source は

\[
 R\in N(W),\quad J=\operatorname{ad}_R|_W,\quad
 T(a)=ta+r_b b+r_c c,\quad T(b,c)=(tI-J)(b,c),\quad3t=\operatorname{tr}J.
\]

determinant (2) の恒等零は、\(tI-J\) が scalar であることと同値。
nondegenerate W では \(R=\alpha H\)、\(J=\operatorname{diag}(-2\alpha,2\alpha)\) なので \(R=0\)。
Borel では \(R=\alpha H+\beta E\)、J の eigenvalues は0,2α、nilpotent off-diagonal は−2βなので、
scalar 条件は同じく \(\alpha=\beta=0\)。
従って全 high-m source は正確に

\[
 R=S=0,\quad T=(r_b b+r_c c)\partial_a,\quad K=0,
 \quad r_b,r_c\in\mathbb C. \tag{3}
\]

旧 m二では除いた R0 全体である。
R0 では旧 row shear \(r\mapsto r-\ell J\) は何も消さないため、全二係数を残さなければならない。

## 4. (3) の actual normalizer と全有限型

**Nondegenerate W。** \(M=bF+cE\) とする。
normal \(\operatorname{diag}(t,t^{-1})\) と base
\(b'=t^{-2}b,c'=t^2c,a'=a\) は U,M を保存し

\[
 (r_b,r_c)\longmapsto(t^2r_b,t^{-2}r_c)
\]

と作用する。normal swap と対応する base swap は両係数を符号付きで交換する。
全非零 row は、\(r_b\ne0\) に交換後 \(t^2=1/r_b\) を選ぶことで

\[
 r=(1,\rho),\qquad\rho=r_b r_c\in\mathbb C
\]

に直る。ρ0 は非零の isotropic axis、ρ≠0 は全 nonisotropic rows。
r0 は別の閉型として残す。
ρ は normalizer invariant なので、全非零 row を無根拠に一点へ置換していない。
平方根は C 上で選べ、逆は t の逆と inverse swap である。

**Borel W。** \(M=bH+cE\) とする。
normal \(g=\exp(sE)\operatorname{diag}(t,t^{-1})\) と base
\(b'=b,c'=-2sb+t^2c,a'=a\) は U,M を保存し

\[
 (r_b,r_c)\longmapsto(r_b+2sr_c/t^2,\ r_c/t^2)
\]

と作用する。
rc≠0 では \(t^2=r_c,s=-r_b/2\) により \(r=(0,1)\) へ可逆に直る。
rc0 では rb は残る invariant であり、\(r=(\rho,0)\)、全ρ∈C を保持する。
後者のρ0は r0 である。
いずれも graded actual base/normal transformation なので、全30/9 higher source は失われない。

従って M(p)0 の全 high-m source は ND(ρ) 一 family、ND-r0 一 fixed 型、
B-rc≠0 一 fixed 型、B-rcline(ρ) 一 family の四つである。

## 5. 消費した parameter を戻した全 source 次元

fixed plane の M(p)非零の全 pair (p,M) は dimension10。
W Borel で一、さらに M0 nilpotent で一を引く。
ND0 は追加 parameter 無し、B0 はτ一を加える。
M(p)0 では pair base は ND8、Borel7、(3) の全 row は二次元。
NDのρ固定 nonzero row は一条件、r0は二条件。
Borel rc≠0 は open、rc0 は一条件、さらにρ固定で一条件である。

| 全 high-m source | affine symbol base | parameter 固定時 | plane6/projectivization−1 を戻した image 上界 |
|---|---:|---:|---|
| NDSS0 | 10 | — | \(15+d\) |
| NDN0 | 9 | — | \(14+d\) |
| BSS0(τ) | 10 | 9 | generic \(15+d\)、fixedτ \(14+d\) |
| BN0(τ) | 9 | 8 | generic \(14+d\)、fixedτ \(13+d\) |
| R0-ND(ρ), row≠0 | 10 | 9（ρ0も同じ） | generic \(15+d\)、fixedρ \(14+d\) |
| R0-ND-r0 | 8 | — | \(13+d\) |
| R0-B-rc≠0 | 9 | — | \(14+d\) |
| R0-B-rcline(ρ) | 8 | 7（ρ0も同じ） | generic \(13+d\)、fixedρ \(12+d\) |

d は全30 weight一 operator の kernel と全9 weight二 operator の kernel の次元和。
任意 P1 で P2 の方程式は同じ linear operator の inhomogeneous system だから、
全 higher Poisson fiber がこの d 以下になる。
高 weight の非線形 Jacobi を解いた、m≥3 が空又は E に帰属したとは主張しない。

NDSS0 は旧 NDSS(s) の s0 に完全一致し、既存 excluded-point certificate の d=8 を復活出来る。
その全像は \(10+8+6-1=23\) で、新しい rank 計算を要しない。
他七型について fixed又はQQ(ρ/τ) の全350行を一回の有限 job で計算し、
selected maximal minor の全 irreducible factor と全 complex roots を exact に保持する。
generic rank だけで exceptional roots を捨てず、固定parameterで一つ下がる base 次元を正確に戻す。
