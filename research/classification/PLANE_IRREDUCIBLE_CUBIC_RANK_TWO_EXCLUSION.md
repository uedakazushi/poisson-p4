# rank二planeの全既約cubic零divisor：node/cuspも含む認識とG排除

2026-11-07、主義務G。smooth cubicに対する既存証明・証明書は変更しない。
本noteは、その独立したnode/cusp拡張である。

**命題。** S=P2⊂P4を実退化因子のSingDに含まれるPoisson invariant planeとし、
π|Sがgeneric rank二、その零divisor Rが既約・被約cubicであるとする。
このとき、rank四πの実Dは既約・被約quinticにならない。
Rのsmoothness、Dの横断node性、blowupのvertical divisor不存在を仮定しない。

## 1. 全既約plane cubicを扱う三つの型

複素数体上の既約plane cubicはsmooth、ordinary node、ordinary cuspのいずれかである。
後二型はそれぞれprojective座標で

`f_node=y²z−x³−x²z`、`f_cusp=x²z−y³`

と取れる。smooth型は既受理の
`PLANE_SMOOTH_CUBIC_CONORMAL_CONNECTION.md` と独立reviewにより処理済みである。
以下ではsingularな二型を全て扱う。

## 2. complementのmeridianは両型でorder三

nodal cubicについては、直接一次資料確認済みの
E. Artal Bartolo–A. Dimca, arXiv:1507.08178v1,
<https://arxiv.org/pdf/1507.08178> §3（printed pp.5–6）を使う。
nodal plane curveのcomplementのπ1が可換であるというFulton–Deligne入力と、
同箇所で明記されたmeridian relation `Σd_i μ_i=0` により、
既約degree三の場合は `π1(P²\R)=Z/3` である。
Fulton–Deligneの原証明を本noteで再証明したとはしていない。

cuspidal cubicでは次の直接van Kampen計算で同じ結論を得る。
U=P2\{x²z−y³=0}とし、`U0=U∩{x≠0}` とする。
`t=y/x, s=z/x, r=s−t³` により `U0≅C×C*` なのでπ1(U0)=Zである。
追加するline L={x=0}からcusp point{y=0}を除いた部分はCであり、y=1 chartでは
`u=x/y, v=z/y`、curve方程式は `u²v−1=0` となる。

U内のL\Rの近傍を
`N={ (u,v): |u|²(1+|v|)<1/2 }`
と取る。Nはuを零へ縮めてCへdeformation retractするのでsimply connected。
`N∩U0` はv上のpunctured disk bundleで、connectedかつπ1=Zである。
v=0におけるnormal meridian u=εexp(iθ)では

`t=1/u, r=(v u²−1)/u³=−u^-3`。

従って `π1(N∩U0)→π1(U0)` はgeneratorを−3倍へ送る。
U=U0∪Nなので、van Kampenにより `π1(U)=Z/3` となる。
追加曲線の近傍による別のgenerator又は関係を無記録で捨ててはいない。

## 3. flat tracefree一位極接続はgenericにlogarithmicとなる

conormal bundleはO(-1)⊗C²で、SingD仮定から `{I,I}⊂I²`。
cotangent algebroid接続は全Sで正則係数を持ち、anchor inverseはf^-1倍なので、
projective tracefree接続は全球

`Ω∈H⁰(P²,Ω¹_P²(R))⊗sl₂`、`dΩ+Ω∧Ω=0`

となる。Rの有限個の特異点を除く一般点で、smooth版と同じ
`Ω=M dx/x+N dy/x`、`[M|R,N|R]=N|R` を得る。
もしN|R≠0ならMの固有値は±1/2で、transverse meridian monodromyの固有値は−1である。
§2のorder三meridianと矛盾する。
従ってΩはRの全一般点でlogarithmicである。

`fΩ` は最初から正則である。
generic log性から `f dΩ` にcodimension一のpoleは無い。
残る候補は有限個の点だけなので、P2のsmoothness/Hartogsにより `f dΩ` も全域で正則となる。
従ってΩの各成分はSaitoの意味のglobal logarithmic一形式である。
特異点でlog性を無条件にnormalizationの正則形式へ置き換えてはいない。

## 4. nodal/cuspidalの全global logarithmic一形式を直接求める

任意のglobal一位極形式αは、homogeneous coneで
`α=i_E i_V(dx∧dy∧dz)/f`
と書ける。Vはlinear vectorfieldで、Eulerの倍を商する。
これは外積Eulerによる `H⁰(Ω¹(3))≅sl₃` という八次元の全sourceである。
αがSaito logarithmicである条件は、正確に `V(f)=λ f`（λ定数）である。
実際、dfとnumeratorのwedgeのmod f classは、V(f)をEuler contractionに掛けたものとなる。
fはprimeなので、fがV(f)を割ることと同値である。

node型では `V=(a x+b y+c z, d x+e y+g z, h x+i y+j z)` とすると、
`V(f_node)=λ f_node` の係数比較から順に
`h=i=b=d=c=g=0, j=e=a, λ=3a` を得る。
従ってVはEulerだけで、`H⁰(Ω¹(logR_node))=0` である。
よってΩ=0。

cusp型では `V(f_cusp)=λ f_cusp` は全off-diagonal係数を零とし、
diagonal係数に `2a+j=3e` を課す。
Eulerを商すると一次元で、代表は `V=(3x,2y,0)`。
対応する全logarithmic一形式は、符号を除いて

`α=(2yz dx−3xz dy+xy dz)/(x²z−y³)`

のscalar倍である。
直接微分すると

`dα=(z dx∧dy−y dx∧dz+x dy∧dz)/(x²z−y³)≠0`。

従ってΩ=α*T（T constant sl2 matrix）と書けるが、
`Ω∧Ω=0` なのでflatnessはdα*T=0、従ってT=0を強制する。
cuspidal型でもΩ=0となる。

## 5. scalar接続の差も消え、全normal-weight証明がそのまま適用する

P_fが与えるcanonical scalar conormal接続との差は、全球閉一形式βのscalar倍である。
βは高々一位極を持つため、上のgeneric log/Hartogs argumentによりSaito log形式となる。
node型では全log形式が零、cusp型では唯一のlog形式αが非閉なので、どちらもβ=0。
従ってconormal接続全体がP_fと一致する。

`PLANE_SMOOTH_CUBIC_RANK_TWO_EXCLUSION.md` のfull30/full9証明は、
fについて**constant Xに対してX(f)=0ならX=0**という性質だけを使う。
既約ternary cubicは二変数だけのhomogeneous cubicではないので、この性質を満たす。
またconstant方向のSchouten及びcurlの全恒等式は任意のternary cubicに成立する。
既存のimmutableな
`CLASSIFICATION_PLANE_SMOOTH_CUBIC_NORMAL_WEIGHTS_20261107.json`
も全十cubic係数をsymbolicに保持した恒等式であり、smoothな一数値例ではない。

従って全nonzero binary quadratic Kについて全tensorは実unipotent shearで
`P_f+K(z,w)∂z∧∂w` へ移り、実Dは `3fK`。
K=0の全sourceはrank二以下である。
これで全既約R型のrank四枝を既約被約quinticから排除した。

## 6. 残す補集合

Rがreducible又はnonreducedなplane、rank零plane、plane以外の曲面は残る。
Rのsingular二型を扱ったことは、D又はそのnormalization/conductor全体を滑らかとしたことではない。
本noteにより既約非正規quintic全体又はG全体が閉じたとは主張しない。
