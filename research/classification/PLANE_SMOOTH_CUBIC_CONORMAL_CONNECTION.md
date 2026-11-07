# rank二plane上のsmooth cubic：conormal接続のlog性を仮定せず証明する

2026-11-07、主義務G。これはplane上のfirst-normal-jetを制限する補題であり、
Poisson全deformation又は未知優越成分をまだ尽くす定理ではない。
blowup逆形式のvertical residueが高位極を持つ可能性は別問題であり、本補題では除外しない。

## 1. flat sl2接続に対する補題

R⊂P2を滑らかな既約三次曲線、U=P2\Rとする。
`Ω∈H⁰(P²,Ω¹_P²(R))⊗sl₂` が
`dΩ+Ω∧Ω=0` を満たすなら、**Ω=0** である。
ここでΩはRに沿い高々一位のmeromorphic matrix一形式という仮定だけであり、
初めからlogarithmicであるとは仮定しない。

Rの局所方程式をx、接方向座標をyとし、
`Ω=M(x,y) dx/x+N(x,y) dy/x`、M,N holomorphic・tracefreeと書く。
flatnessのx^-2係数は

`[M(0,y),N(0,y)]=N(0,y)`。                              (1)

もしΩがlogarithmicでなければ、Rの一般点でN(0,y)≠0である。
式(1)から `tr(N²)=0`、従ってNは非零nilpotentである。
NをE12へ共役すると、(1)はMの下三角成分が零、二対角成分の差が一であることを強制する。
traceM=0なのでMの固有値は `1/2,-1/2` である。

その点を通るtransverse diskでyを固定すると、接続d+Ωは
`d+M(x,y0) dx/x` という通常のregular singular接続である。
regular singular一次系のresidue–monodromy spectrumにより、meridian monodromyの固有値は
`exp(±2πi/2)=−1` である。resonanceがあっても固有値のこの結論は変わらない。

一方d+ΩはUの**自明rank二bundle上の全球flat SL2接続**である。
smooth degree三plane curveのcomplementは `π1(U)=Z/3` なので、
meridian monodromyの三乗はidentityでなければならない。
固有値−1はこれに反する。
従ってN|R=0、すなわちΩはlogarithmicである。

residue完全列
`0→Ω¹_P²→Ω¹_P²(logR)→O_R→0`
では、`H⁰(O_R)=C→H¹(Ω¹_P²)=C` が1をc1(R)=3Hへ送り、単射である。
`H⁰(Ω¹_P²)=0` と合わせて `H⁰(Ω¹_P²(logR))=0` を得る。
よってΩ=0となる。

## 2. 使ったtopology/ODE入力の射程

smooth plane curve degree dのcomplementのπ1がZ/dとなる標準定理を使う。
この特殊場合は次のようにも導ける。
R={f3=0}のaffine coneは原点だけをcritical pointに持ち、
homogeneous Milnor fiber `{f3=1}⊂C³` はisolated hypersurface singularityのMilnor fiberと同じhomotopy typeを持つ。
J. Milnor, *Singular points of complex hypersurfaces*,
Annals of Mathematics Studies61, Princeton University Press(1968), Theorem6.5により、
これは二次元sphereのbouquet、従ってsimply connectedである。
自然写像 `{f3=1}→P²\R` は自由なμ3作用のquotientなので、π1(U)=μ3となる。
これは本noteにおける明記した標準topology入力であり、有限行列計算で検証した事実ではない。

直接読んだ補助照合はE. Artal Bartolo–A. Dimca,
*On fundamental groups of plane curve complements*, arXiv:1507.08178v1,
<https://arxiv.org/pdf/1507.08178> のIntroductionと§3、Proposition3.1の証明である。
そこには上記Milnor fiberのcyclic d-fold quotient、nodal curveのcomplementの可換性、
及び可換な場合のmeridian relation `Σd_i μ_i=0` が明記されており、
滑らかな既約三次の場合にZ/3と一致する。
Milnor原著の全証明を本作業で再証明したとはしていない。

もう一つの標準入力は、disc上のholomorphic residue Mを持つregular singular接続で
monodromyの固有値がexp(−2πi eigen(M))となる局所ODE定理である。
非logarithmicな二変数接続全体のregular singular性は用いず、transverse diskだけに適用する。

## 3. Poisson conormal接続へ適用するための一位極の証明

S=P2⊂P4をsmooth Poisson invariant planeとし、π|Sがgeneric rank二とする。
さらに**Sが実退化因子DのSingDに含まれる**と仮定する。
本研究の非正規曲面成分はこの仮定を満たす。
その零divisor Rをsmooth irreducible cubicと仮定する。
`N*=I_S/I_S²=O_P²(-1)⊗C²`。
Sのgeneric rank二点ではWeinstein分解によりπは接方向symplectic二次元とnormal Poisson二次元の積となる。
実Dの一次normal symbolはnormal-normal bracketの一次symbolそのものなので、
`D∈I_S²` は ` {I_S,I_S}⊂I_S² ` を強制する。
このbundle bracketのgenericな消滅は被約Sの全域へ延長する。
従って

`∇_(df)(j)={f_tilde,j_tilde} mod I_S²`

はN*上のwell-definedなcotangent Lie algebroid接続である。
全coefficientはS上で正則で、Jacobiによりflatである。
Uではanchor π_S^sharpが可逆なので、ordinary flat connectionへ移せる。
πSはRで一位に零となるため、そのinverseは一位極だけを持ち、
この接続の一形式もRに沿い高々一位の極を持つ。

P(N*)=P2×P1には正準の自明projective frameがある。
N*の二frame間のtransitionは同じO(-1)のscalarなので、
tracefree partは全球

`Ω∈H⁰(P²,Ω¹_P²(R))⊗sl₂`

にglueする。projective flatnessとtracefree性から `dΩ+Ω∧Ω=0` である。
従って§1を適用してΩ=0を得る。

**Poisson系への結論。** この全plane枝のconormal接続はprojectively scalarである。
constant normal directionは全てprojectively horizontalとなる。
log性は仮定ではなく、上のglobal monodromyにより証明した。
X、導手、Dの横断ordinary-double又はvertical divisorに関する仮定を、このfirst-jet結論では使わない。
SがSingDに入るという仮定は必要であり、単にPoisson invariantでrank二というだけから
`{I_S,I_S}⊂I_S²` を推論してはいない。

## 4. 残る仕事と誤用を避ける境界

この補題はhigher normal-weight項の全Jacobi/Pfaffianをまだ解いていないので、G排除を意味しない。
normal-normal quadraticが二相異なる根を持つ場合だけ、GL2 normal変更でa*z*w型へ直せる。
平方又は零のnormal quadraticを、generic doubleという追加仮定なしに捨ててはならない。

Rがsingular又はreducibleならπ1とresidue完全列が変わり、この証明を適用していない。
πSの零divisorが非被約ならinverse接続は高位極を持ち得るため、同様に補集合である。
blowup inverse二形式のE-residue βがvertical曲線に沿い高位極を持つ場合も、
ここで証明したconormal接続の一位極と混同していない。
