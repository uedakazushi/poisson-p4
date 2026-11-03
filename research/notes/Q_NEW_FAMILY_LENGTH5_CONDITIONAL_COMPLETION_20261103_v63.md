# 条件付き補題：generic 二次 rank 上界から長さ5と完成環を得る

2026-11-03。これは**条件付きの証明原稿**である。
generic な二法線 Hessian の cokernel rank が1以下であるという追加仮定を、この原稿では検証していない。
その exact 証明が独立受理されるまで、以下の自由性・長さ5・完成環の結論を無条件の成果へ昇格しない。

## 入力と未検証の仮定

n111における全 canonical Poisson 完備局所環を A とする。
既証の actual 三変数族・群・scalar の像の26方向から、実 smooth formal branch と、
R=C[[s1,...,s26]]→A→R という恒等合成を得る。全350式の実二法線閉 fiber の独立監査により

    A/m_R A ≅ C[u,v]/(uv,u³,v³)

で、長さ5、Hilbert function (1,2,2)、二次 relation の discriminant は非零である。
完備 Nakayama から A は R 上有限で、R^5→A という全射を持つ。
ここで26方向はすべて群方向ではなく、この時点で A=R⊗B や平坦性は仮定していない。

新しい実成分の次元は canonical affine26・projective25と確定している。
generic Jacobi rank98の独立証明により、この成分に沿う generic 横断 Artin 環の埋め込み次元は2である。

追加で必要な未検証仮定(H)は次である。

    同じ generic 横断二法線における二次 initial relation 空間の rank が1以下。

普遍 polynomial tangent 二本を使う場合には、全 Jacobi 像との augmented rank が99以下であることが
その証明になり得る。ただし actual family26方向との quotient が正しいことも確認する。
この原稿は(H)を数値標本から推論しない。

## 有限 rank の比較は、平坦性を仮定しない

K=Frac(R)と置く。augmentation A→R が選ぶ generic branch は residue field K を持つ。
その generic Artin 局所環の極大 ideal を n とすると、埋め込み次元2と(H)により

    dim K(n/n²)=2,   dim K(n²/n³)≥3−1=2

である。したがってその長さは少なくとも1+2+2=5となる。
この branch の長さだけで dim_K(A⊗K)≥5を得る。
一方 R^5→A の全射から dim_K(A⊗K)≤5であり、両者は5に等しい。

全射の kernel M は R^5 の部分加群で、Rが整域なので torsion-free である。
M⊗K=0から M=0となり、Aは**この段階で初めて** R上自由 rank5と結論できる。
最初からAのtorsion-freenessやflatnessを仮定する必要はない。
この比較により他の generic branch の rank も残らず、augmentation が選ぶ局所因子が全A⊗Kを使う。

単に閉 fiber の長さが5であることだけから kernel Mを零にしてはいけない。
kernel の閉 fiber への写像には Tor の問題があり得る。
上の generic rank 下界5こそが、その誤った推論を避ける追加の根拠である。

## nilpotent ideal と乗法

N=ker(A→R)と置く。augmentation が split なので A=R⊕N であり、Nは自由 rank4である。
generic Artin 環は長さ5、埋め込み次元2、二次部分の次元が少なくとも2なので、
その Hilbert function は厳密に(1,2,2)、N_K³=0となる。
AのR自由性から N³=0を全完成環へ戻せる。

N/N²は閉 fiber で二生成、generic rank2である。
完備 Nakayama による R²→N/N² の全射の kernel も自由加群の rank0部分加群なので零であり、
N/N²は自由 rank2となる。この商への exact sequence は R加群として split し、N²も自由 rank2である。
u,vを N/N²の基底の lift とすると、N³=0から N²はu²,uv,v²で生成される。
従って Sym²(R²)→N² の kernel は自由 rank1で、ある二次形式qを生成元に取れる。
閉 fiber の非退化性からqの discriminant は Rの単元である。

Rは複素数上の完備局所環なので、Henselの平方根と可逆な二変数線形変換によりqをuvへ移せる。
この時

    R[u,v]/(uv,u³,v³) → A

は全射であり、両辺ともR自由 rank5、閉 fiber でも同型なので同型となる。
したがって(H)の受理後には、actual canonical 完成環を

    C[[s1,...,s26,u,v]]/(uv,u³,v³)

と同定できる。非零 cone の自由 scalar 方向を一つ除き、projective 完成環は

    C[[s1,...,s25,u,v]]/(uv,u³,v³)

となる。この完成環には極小素が一つで、被約環は滑らか、nilradicalの三乗は零である。
Artin因子の socle 次元は2なので、全完成環は Cohen–Macaulay・type2で、Gorensteinではない。
これらもすべて(H)を必要とする条件付き結論である。

## 二次 obstruction の quotient に関する注意

追加 tangent Tが actual family26に横断する universal section なら、
J(P(t))*T(t)=0をパラメータ方向で微分して

    Bilinear(P'(t),T(t)) + J(P(t))*T'(t)=0

を得る。従って family方向との mixed二次項はJacobi像で消える。
群・scalar方向では同じことが有限作用の微分と同変性から従い、family方向自身の二次は
全次数実族の二階微分で積分される。
これにより universal 二法線の Hessian cokernel は26方向を割った intrinsic 二法線 obstruction を測る。
一点だけで選んだ arbitrary kernel vector の Hessian に、この降下を無条件で適用してはいけない。

## 現在の状態

actual 成分25、generic Jacobi rank98、二つの universal polynomial tangent、
および n111 の実閉 Artin fiber長さ5は、それぞれ別の凍結証明・独立 review にある。
本原稿は有限-rank補題の適用と完成環の代数的帰結を監査したもので、
(H)の証明やその受理を代行しない。
