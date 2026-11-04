# 実 CA-critical 接触模型、比較 kernel rank4 と残る深さ判定

2026-11-04。対象 D は **実際の scheme 交叉 CA∩CB** であり、全 Poisson
Jacobi scheme そのものではない。V1 は以前から保存されている CA-critical V
稠密軌道上の同じ有理点である。新しい係数計算は実 source 三制限の3300 slot
だけで完了し、source と純代数の二つの独立監査がともに PASS した。
本稿はその採用 corollary であり、B18 原点の深さ自体はまだ決めない。

## 1. V1 における全交叉と実比較写像

G=k[[g1,...,g20]] とする。全150出力を残した実 source 計算、N に適合した
fold、CB の単位 even29 射影、および全 ambient graph の接触係数から、

    D_V = G[[l1,...,l7,v,x]] / (x*l1,...,x*l7,x^2,x*v^2)

を得る。ここでは七つの l と v,x は transverse 九座標である。実際の二成分の
normalizer を同じ ambient 上で tensor した比較対象は

    E_V = G[[l1,...,l7,v,t]] / (l1*t,...,l7*t,t^2-v*t,t^3),
    D_V -> E_V : x -> t^2,  l_i -> l_i,  v -> v.

この写像の kernel と cokernel はそれぞれ

    K_V = G*(x*v),            C_V = G*class(t).

従って両者の一般 transverse 長さは正確に1である。kernel の長さを二つの
nil-module の長さの差だけから推測したのではなく、上の同じ実写像で計算した。

証明の接触部分を要約する。実 CA の Ac8 座標を
(u8,x=p^2,Y_i=u_i*p,Z=p^3) とし、p は実 N divisor を定める proper-source
第三 binary 係数そのものと取る。これは旧 source 座標22であり、fold 微分は1。
選んだ u8 は実成分 A に属する八つの係数函数である。unit fold rescaling は
conductor (u8,p^2)B を保つので、この p を使っても実 Ac8 ring は変わらない。
CB の同じ even 座標への射影行列式は128/20599365234375で非零であり、
N は CB 上で x=0、したがって Y_i=x*b_i、Z=x*c と書ける。
実 row60 の b_i(0)=20/9 は単元である。

v=u_i/b_i と置いて七つの l_j=u_j-v*b_j を使うと、全81 intrinsic 関係の
x を除いた ideal は厳密に (l7,x-v^2,v*(v-c)) に等しい。追加の全 ambient
graph は x*G_j の形で保持する。八つの単位式を解いた後、実 source 計算は
c の v 一次係数と全 G_j の v 一次係数が零であることを証明した。
従って ctilde と Gtilde_j は全次数で (v^2) に入り、v*(v-ctilde) は v^2
の単元倍となる。高次 graph はこの長さ2商を変えられない。これが D_V の
全次数 ideal の証明である。

さらに t^2=x、u_i*t=x*b_i、x*t=x*c を付加する実 normalizer base change では
v*t=t^2、l_i*t=0、t^3=0 を得て、逆向きの全式もこれから従う。全 graph は
そこで零となる。これが E_V と上の実比較写像を与える。

一点での係数零を、無関係な20変数函数の零と扱わない。V1 の群軌道微分は20
であり、群の20次元 section と canonical ambient の横断補空間への作用は
形式的 IFT 同型を与える。CA、CB、N とその proper normalizer 比較写像を
同時に引き戻すと共通 G20 product になる。従ってこの横断係数の判定は V1
稠密軌道上の全 diagram に一様に適用される。

Nil(D_V) は G-free2、basis x,xv、square-zero である。被約商は正則28次元、
D_V の embedding dimension は29、depth は20である。projective 版では
共通 scalar を除き dim27、edim28、depth19となる。これは V1 の局所結論であり、
B18 原点 P0 の深さを決めるものではない。新しい全域 associated closure や
Poisson 成分の追加数は本稿では主張しない。

## 2. 同じ零weight比較への厳密な transport

B18 共通 G19 を除いた相対 whole graded ring を D^gr、その同じ二成分の
normalizer tensor を E^gr とし、K^gr=ker(D^gr->E^gr) とする。実共通底は

    R10=R0[[w]],  A_i=R10+I_i*tau_i subset B_i=R10+R10*tau_i,
    I_i=I_i0+w*J_i0,            deg(tau_i)=0.

これは共通 physical 底での同じ実埋込みであり、別々の抽象 normal form を
同一視していない。f in I_A、g in I_B の乗法は、両 normalizer から両実成分へ
戻る ambient-linear 写像を与える。その二つを tensor すれば比較写像との
両合成は fg 倍なので、I_A*I_B は K^gr と比較 cokernel を annihilate する。
flat tensor や Tor 消去を使わない。

I_A0 の一般点では I_B0 は単元である。実 I_A0+I_B0 商が有限長だからである。
また I_A0+J_A0 は実 F9 maximal ideal なので J_A0 もその一般点で単元。
従って w belongs to I_A であり、この一般比較 kernel は w に殺される。
反対側の二枝にも同じ議論を使える。

ゆえに各一般 kernel は実 w-torsion T(w) に入る。受理済み T1=T2=T3=0 と
degree>=4 の安定化は、T(w)_n=0 for ALL n>=1 を意味する。非負 grading の
exact degree0 部分を取ると

    (D^gr)_0=D0,  (E^gr)_0=E4,
    K=ker(D0->E4)

となり、上の一般 kernel の全長がこの同じ K に残る。w に殺されるだけで
degree0と推論したのではなく、正 degree の全 torsion 消滅を用いている。
また K と ker(D0->D4) を同一視していない。共通形式的 G19 product は flat
なので、全比較 diagram と一般長を保つ。

## 3. K の正確な DVR rank

受理済み full conductor/source inventory により K の一般支持として調べる
枝は CA-critical、CA-selfpair Pa、CB-critical、CB-selfpair Qb の四本で尽くされる。
R の一般点では両 conductor は単元なので比較は同型である。実 Noether
parameter ell は四本すべての reduced branch 上で order1であり、枝の一般
定数体も k のままである。DVR を D=k[[ell]] と書くと、一般長の寄与は

| 実 source 枝 | K の一般長 |
| --- | ---: |
| CA-critical V | 1 |
| CA-selfpair Pa | 1 |
| CB-critical V | 0 |
| CB-selfpair Qb | 2 |

Pa は受理済み full D=K[[u8,f]]/(f^2,f*u8) と実 normalizer 比較から1。
Qb は full D=B/(f*J_Z)、J_Z=(z7,f^2,rho^2) と実比較から
((f) intersect J_Z)/(f*J_Z)=span(f^2,f^2*rho) で2である。
CB-critical は実 CA の外にある。Pa/Qb の scheme triple-meeting はそれぞれ
長さ1/2であり、被約支持だけの議論には置き換えていない。したがって

    rank_D K = 1+1+0+2 = 4.

これは K の閉点 torsion の長さや、その消滅の結論ではない。

## 4. 未計算の有限 T と原点の深さを結ぶ数式

同じ実比較の完全 complex を R=R0 上で

    R^752 --X--> R^40 --V--> R^3

とし、T=ker(V mod ell)/im(X mod ell) とする。X の完全性には受理済みの
actual chain lift と完全 normalizer syzygy を用いる。T を単純に20次元底上の
未修正 kernel と定義してはいけない。受理済みの正確な特殊化列は

    0 -> K/ell*K -> T -> C[ell] -> 0.

C=coker(D0->E4) は finite D-module、mu_D C=8、rank_D C=4 なので
dim_k C[ell]=4。K も conductor 積に殺され finite D-module であり、T は
有限長である。従って未計算の L=length_k T に対し

    mu_D K = L-4,
    mu_D(tors_D K) = mu_D K-rank_D K = L-8,
    L >= 8.

ここで torsion の最小生成数と torsion の長さ・指数を区別する。
DVR 上で K=D^4 direct_sum tors_D K なので、L=8 は K が D-free4であることと
同値である。L>8なら torsion は非零だが、その長さはまだ分からない。

正規化比較 E4 は受理済み実 source CI 模型により H0_m(E4)=0、depth1。
従って D0->E4 の像にも H0_m はなく、H0_m(D0)=H0_m(K) である。
K の conductor 支持上では ell が Noether parameter なので、後者は
tors_D K に等しい。さらに全 positive D_n は受理済み T1=T2=T3=0 により
深さ1の安定 tail D4 に一致する。従って whole graded D^gr の m-socle は
degree0 の H0_m(D0) に正確に一致する。

共通 G19 product と元の深さ上下界 {19,20} を戻せば、B18 原点において

    depth_cone(actual D at P0)=20  iff  L=8,
    depth_cone(actual D at P0)=19  iff  L>8.

projective 深さはそれぞれ19/18である。これは**完全に定義された有限 T の
未計算長 L による判定定理**である。本稿は L、閉点 torsion、原点の深さや
全 nil 構造を計算していない。20次元底への安易な kernel base change、348底
だけからの有限 jet 推測、相対 CA/CB 全域分類の主張も行わない。

必須 peer は CONTACT (actual source)、CONTACT_ALGEBRA (全接触代数・zero-weight
transport)、COMPARISON_C_GENERIC_DVR_RANK、SPECIALIZED_COMPARISON_T_PREFLIGHT、
T1_T2_QUADRATIC_PHI、T3_ZERO_STRUCTURAL と proper core CI の各 REVIEW。
全ての exact path と SHA256 を同名 report に束ねる。
