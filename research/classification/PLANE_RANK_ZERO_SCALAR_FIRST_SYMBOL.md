# rank 零 plane・traceless co-Higgs 零：全 scalar first symbol

2026-11-07。主義務 G、producer `plane_collision_review`。
実既約被約 quintic D を持つ rank 四 Poisson tensor が plane P 上零となるとする。
新しい全重複度 conormal 可換性の下で、traceless co-Higgs が **恒等的に零** となる全 first symbol を扱う。
ここで「scalar」は非零 traceless part と trace field が pointwise dependent である場合を含めない。
その global-V dependent source は別の未完了義務である。

以下は pure first-jet 幾何・全 actual normalization と dimension count。
finite higher-operator bound は別に実行し、未検査モデルを G 完了とはしない。
T0K非零は別作者の手証明、T0K0は別残余として明示する。

## 1. 全35源からの正確な scalar 11源

`REDUCED_DEGENERACY_ZERO_SURFACE_ABELIAN_CONORMAL.md` と別作者 review より、
rank 零 plane の conormal Lie bracket は被約実 D の全可能重複度で零。
従って既受理の full35/30/9 normal-weight source を使える。
35の weight零 source において traceless co-Higgs の全三 entry を零とすれば、残る全源は

\[
 \Pi_0=T\wedge E_n+K(z,w)\partial_z\wedge\partial_w,
 \quad E_n=z\partial_z+w\partial_w,\quad
 T\in\mathfrak{sl}_3,\quad K=q_S,\ S\in\mathfrak{sl}_2.
\]

\(q_S=\det((z,w),S(z,w))\) とし、
\(S=\left(\begin{smallmatrix}s&t\\u&-s\end{smallmatrix}\right)\) なら
\(q_S=uz^2-2szw-tw^2\)。全 binary quadratic と一対一である。
T の trace-free lift は base Euler の gauge を除く唯一の代表。
従って全 scalar first-symbol affine space は正確に \(8+3=11\) 次元であり、
特定の diagonal T から全源を推論したものではない。

全 T,S について Jacobi は零である。
T は normal variables に依らず、\([T,E_n]=0\)、
\([E_n,K\partial_z\wedge\partial_w]=0\)（normal weight零）なので全 Schouten cross term が零。
さらに \(\Pi_0^2=0\)。T∧En と vertical bivector の wedge は normal dimension 二により零だからである。
従って leading actual D2 は常に零で、scalar first symbol は m二には存在せず、全 high-m へ送られる。
high-m の空性をこの事実だけから推論しない。

base g、normal h、tensor scalar σ の作用は正確に

\[
 (T,S)\longmapsto(\sigma gTg^{-1},\ \sigma hSh^{-1}). \tag{1}
\]

normal determinant factor は q の変換と \(\partial_z\wedge\partial_w\) の変換で相殺する。
全変換は grading を保つ actual PGL transformation と非零 scalar であり、
全30 weight一・全9 weight二 unknown を保ち、その逆で元の全源へ戻る。

## 2. T≠0 の全五 Jordan 層と全 parameter 次元

\(D=\operatorname{diag}(1,1,-2)\)、\(N_2=E_{01}\)、\(N_3=E_{01}+E_{12}\) とする。
trace零3×3行列の非零 Jordan 型は次の五つだけ。

| T の全型 | 代表 | affine T-stratum 次元 |
|---|---|---:|
| regular semisimple | diag(λ,μ,−λ−μ)、三固有値相異 | 8 |
| repeated semisimple | λD、λ≠0 | 5 |
| mixed | λD+N2、λ≠0 | 7 |
| regular nilpotent | N3 | 6 |
| rank一 nilpotent | N2 | 4 |

三固有値相異なら diagonalizable。
二固有値なら repeated block が diagonal 又は size二 Jordan の二場合。
一固有値なら trace零により全固有値零で、非零 partition は3又は2+1。
従って T0を除き全てを尽くす。
次元は GL3 orbit の \(9-\dim Z_{GL_3}(T)\) に eigenparameter 数を足したもの。
centralizer 次元は順に3,5,3,3,5である。

S の全三型は semisimple \(\kappa H,\kappa\ne0\)、nonzero nilpotent E、zero。
\(H=\operatorname{diag}(1,-1)\)、\(E=E_{01}\)。
その affine source dimensions は順に3,2,0。
対応する fixed K は \(-2zw,-w^2,0\) とする。
零・nilpotent・semisimple を一つの generic orbit へ混ぜていない。

## 3. 全 tensor/residual scaling とその逆

S semisimple なら normal conjugation で κHとし、tensor \(\sigma=1/\kappa\) で S=Hにする。
このとき T の eigenvalues は κで割られ、全 residual eigenparameter を保持する。
κの符号・Jordan basis の finite choice は dimension を増やさない。

S nilpotent なら nonzero scalar multiple を normal conjugation で Eへ直せる。
さらに tensor σを掛けた後も、normal diagonal \(h=\operatorname{diag}(t,t^{-1})\)、
\(t^2=1/\sigma\) を選べば \(h(\sigma E)h^{-1}=E\)。
S0は任意の非零 σで0のまま。
従って K0と K-nil の両型では T の非零 eigenvalue を一へ正規化出来る。

nilpotent T は tensor σでの係数を actual base diagonal で復元する。
N3の場合 \(g=\operatorname{diag}(\sigma^{-2},\sigma^{-1},1)\)、
N2の場合 \(g=\operatorname{diag}(\sigma^{-1},1,1)\) が
\(g(\sigma N_i)g^{-1}=N_i\) を満たす。
mixed の N2係数も同じ後者の gで一へ戻せる。
全逆は g,h,σの逆であり、平方根が必要な場合は C上で選ぶ。
unknown higher-weight coefficient を制約する操作は無い。

## 4. T≠0 の全十五 model source

以上により、次の10固定・4一parameter・1二parameterが全 T≠0を覆う。

**固定十型。** N3又はN2と S0,E,Hの全六組、
T=D又はD+N2と S0,Eの全四組。

**一parameter四型。**

* T=λD、S=H、全λ≠0。
* T=λD+N2、S=H、全λ≠0。
* T=diag(1,t,−1−t)、S0。
* 同じ T、S=E。

後二つの域は \(t\notin\{1,-2,-1/2\}\)。
Tが非零 regular semisimple なら有限の eigenvalue permutation により一つの非零 eigenvalue を選べ、
tensor scalingで一へ正規化するので、この一parameter source は全源を覆う。
t0,−1の零固有値は許し、除外三値は repeated semisimple sourceへ送る。
λ0の補集合は T0又はN2として既に別に保持した。

**二parameter一型。** T=diag(λ,μ,−λ−μ)、S=H、
\((\lambda-\mu)(2\lambda+\mu)(\lambda+2\mu)\ne0\)。
ここでは λ0やμ0を捨てず、三固有値相異だけを課す。
一つの点の rank を全二parameterへ拡張してはいない。

## 5. affine source を戻す弱い G 目標

fixed plane の全 scalar symbol dimension は \(d_T+d_S\)。
全30/9 operator による higher-fiber upper を d とすると、plane choice6・tensor projectivization−1を戻した
全像上界は \(d_T+d_S+5+d\) である。
normalization で消費した eigenvalue/scale はこの affine source 次元へ既に戻している。

| fixed model | affine symbol source | G offset | G≤24に十分な d |
|---|---:|---:|---:|
| N3,S0 | 6 | 11 | 13 |
| N3,E | 8 | 13 | 11 |
| N3,H | 9 | 14 | 10 |
| N2,S0 | 4 | 9 | 15 |
| N2,E | 6 | 11 | 13 |
| N2,H | 7 | 12 | 12 |
| D,S0 | 5 | 10 | 14 |
| D,E | 7 | 12 | 12 |
| D+N2,S0 | 7 | 12 | 12 |
| D+N2,E | 9 | 14 | 10 |

一parameter源の全 affine symbol 次元は順に8,10,8,10、G offsetは13,15,13,15。
fixed parameter の有限 wall では一つ下がる。
二parameter regular-ss/H源の全 symbol dimension は11、offset16。
generic、真のcurve、finitepointではそれぞれ d≤8,9,10が十分。
全 minor factor/root と全補集合を保持するまでは、その源の完了とはしない。

## 6. T0と他 source の境界

T0K非零では generic_exhaustion の全 source 手証明が別に進行している。
Kzwのactual shear核とKz2の残るdecomposable fieldを使うもので、本ノートでは重複計算しない。
T0K0は全weight零symbolが零となる別残余であり、上の表へ無証明で入れない。
traceless co-Higgs 非零で trace と pointwise dependent な global-V源も別に残す。
今回の scalar 11源とこの dependent源を同一視しない。
本ノートは全first-jet還元の証明で、有限higher bounds前に全Gを閉じたとは主張しない。
