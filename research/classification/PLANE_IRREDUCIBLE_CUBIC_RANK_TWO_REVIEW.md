# 全既約cubic零divisorを持つrank二plane：node/cusp拡張の独立監査

2026年11月7日、主義務 `E_reduced_nonnormal`（Gの残余排除）。
**判定：記載された仮定でPASS。**
Poisson不変plane S⊂SingDについて、π|Sが一般rank二で零divisor Rが既約・被約cubicなら、
rank四πの実Dは既約・被約quinticにならない。
Rのsmoothness、Dの横断ordinary-double、vertical divisor不存在は追加しない。

対象producer noteは `PLANE_IRREDUCIBLE_CUBIC_RANK_TWO_EXCLUSION.md`。
独立有限証明書は `CLASSIFICATION_PLANE_IRREDUCIBLE_CUBIC_RANK_TWO_REVIEW_20261107.json`、
コードは `classification_plane_irreducible_cubic_rank_two_review.py`、
契約は `jobs/plane_irreducible_cubic_rank_two_review.json`。
Fractionのみの全係数計算と、以下の幾何的な別読解を区別する。
元のsmooth接続review、全f十係数を保持したfull30/full9reviewは不変の既済入力とする。

## 1. 残る全型とnodal complementの外部入力

integral plane cubicの算術種数は一で、singularならnormalizationとの差δの和も一。
従って一つのnode又はcuspだけを持ち、複素PGL座標で

\[
 f_{\mathrm n}=y^2z-x^3-x^2z,\qquad
 f_{\mathrm c}=x^2z-y^3
\]

に取れる。smooth型は既済証明を用いる。

nodal型については、E. Artal Bartolo–A. Dimca,
[*On fundamental groups of plane curve complements*, arXiv:1507.08178v1](https://arxiv.org/pdf/1507.08178v1)
の§3、printed p.5（PDF index4）を、このreviewで直接読んだ。
そこではnodal curveのcomplement群の可換性をFulton–Deligneの結果として明記し、
Proposition3.1の証明で、可換群は各既約成分のmeridianで生成され次数を係数とする一つの関係を持つとしている。
従って本件の既約degree三ではπ₁(P²∖R)=Z/3で、smooth locusのmeridianがgeneratorとなる。
references[5],[10]もprinted p.6で照合した。
Fulton又はDeligneの原証明を今回直接読解したとは主張しない。

## 2. cuspのtwo-open計算を直接点検

U=P²∖{x²z−y³=0}、U₀=U∩{x≠0}とする。
t=y/x、r=z/x−t³によってU₀=C×C*で、rの正のloopがRのmeridianである。
U∖U₀はx=0上のcusp以外の点、すなわちy≠0の一つのaffine lineだけである。
y=1 chartでu=x/y、v=z/yと書き、

\[
 N=\{(u,v):|u|^2(1+|v|)<1/2\}
\]

と置くと、|u²v|<1/2だからNはUに入り、U₀∪N=Uである。
uを零へ縮めるdeformationはN内に留まり、NはCへretractしてsimply connectedとなる。
N∩U₀はu≠0で、v∈C上のpunctured disk bundle。
正の連続radius (2(1+|v|))⁻¹ᐟ²でtrivializeでき、connectedでπ₁=Zである。
v=0、u=εexp(iθ)のgeneratorに沿って

\[
 r=(vu^2-1)/u^3=-u^{-3}
\]

なので、intersectionからU₀への写像はgeneratorを−3倍へ送る。
van Kampenによりπ₁(U)=Z/3であり、U₀のmeridianの像もorder三。
追加したline全体とそのintersectionを含めて計算しており、cusp近傍を無記録で捨てていない。

## 3. tracefree接続のgeneric log性とsingular点での延長

既済smooth reviewのbundle・flatnessの議論はRのsmoothnessを用いない。
S⊂SingDとgeneric rank二から{I,I}⊂I²を得、conormal接続は全Sで定義される。
N*=O(−1)⊗C²なのでtracefree部分は全球

\[
 \Omega\in H^0(P^2,\Omega^1(R))\otimes\mathfrak{sl}_2,
 \qquad d\Omega+\Omega\wedge\Omega=0
\]

としてglueする。特にU上のd+Ωは自明rank二bundleのhonestなflat SL₂接続である。
projective monodromyを未証明のままSL₂へliftしてはいない。

Rのsmoothな一般点でlocal equation tと接座標sを使い、
Ω=M dt/t+N ds/tと書く。
flatnessのt⁻²係数は[M₀,N₀]=N₀。
nonlogならN₀≠0、tracefree性からM₀の固有値は±1/2となる。
sを固定したtransverse circleを零へ縮めると、monodromyの固有値はともに−1となる。
円の半径の変更はradial平行移動による共役なので、この固有値計算に非共鳴仮定は不要である。
§1–2のorder三meridianに対するSL₂ monodromyではこれは不可能である。
従ってΩはRのsmooth locusでgenericにlogarithmicとなる。

global一位極は最初からfΩが正則という意味である。
generic log性からf dΩにも全codimension一でpoleが無い。
残る候補はRの有限特異点だけだが、P²はsmoothでΩ²はlocally freeなのでHartogsにより延長する。
よってΩの各成分はSaitoの意味でglobal logarithmic、すなわちfαとf dαが正則である。
特異点上でlog形式をnormalization上の正則形式へ置き換えてはいない。

## 4. 全global log一形式の完全八次元sourceを独立計算

Euler列からh⁰(Ω¹(3))=8。
全形式はlinear field VをEuler倍で商した

\[
 \alpha=i_Ei_V(dx\wedge dy\wedge dz)/f
\]

で尽くされる。dfとのwedgeはmod fでV(f)のEuler contraction倍となる。
fはprimeかつdegree三なので、Saito log条件は正確にV(f)=λf（λ定数）である。

独立コードはVの九係数とλの全十係数から、三次の全十単項式に対する行列を構成した。
全行列・完全kernelはreportに保存した。

| cubic | 全十列行列rank | kernel | Euler商後のlog形式 |
|---|---:|---|---:|
| node | 9 | Eulerのみ | 0 |
| cusp | 8 | Eulerと(3x,2y,0) | 1 |

従ってnodeではΩ=0である。
cuspの唯一のlog形式はscalar倍と符号を除いて

\[
 \alpha=\frac{2yz\,dx-3xz\,dy+xy\,dz}{x^2z-y^3}
\]

であり、全分母を払って直接微分すると

\[
 d\alpha=\frac{z\,dx\wedge dy-y\,dx\wedge dz+x\,dy\wedge dz}{x^2z-y^3}\ne0.
\]

この全三成分の恒等式もFractionで独立確認した。
Ω=αT、T∈sl₂ constantと書け、Ω∧Ω=0だからflatnessはdαT=0を強制してT=0となる。
cuspでもtracefree接続は零である。

## 5. canonical scalar接続との差と全normal-weight源への帰結

P_fのscalar conormal接続との差は、global一位極の閉一形式β·idである。
閉性からRのsmooth locusでgeneric logとなり、同じHartogs議論でSaito logとなる。
nodeでは全log形式零、cuspでは唯一のαが非閉なので、いずれもβ=0である。
従ってconormal接続全体がcanonical接続に等しい。
元のO(−1)のscalar接続自体が零と結論したのではない。

既済full30/full9証明がfに使った性質は、constant Xに対してXf=0⇒X=0だけである。
任意のintegral ternary cubicが二変数だけなら、複素数体上のbinary cubicとして一次式に分解して矛盾する。
従ってこのinjectivityはnode/cuspにも成立する。
両normal formについて全三partial derivativeの独立性も正確計算で照合したが、
全型への含意はこのirreducibilityの議論による。
既済のuniversal f十係数のSchouten・Euler・有限shear恒等式は特殊化せずそのまま適用できる。

全nonzero Kについて実shear後はP_f+K∂z∧∂w、実D=3fKとなり、既約被約quinticではない。
K=0の全源はbase三方向に入りrank二以下である。
したがって全既約・被約cubic零divisorを持つrank二plane枝の排除が閉じる。
reducible又はnonreducedなcubic、rank零plane、非plane曲面を含めておらず、G全体の空性も主張しない。
