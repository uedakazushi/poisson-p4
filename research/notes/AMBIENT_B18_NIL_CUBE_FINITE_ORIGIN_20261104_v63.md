# B18 実交叉の nil 三乗は零weight原点の有限長障害である

2026-11-04。対象は実交叉 D=CA∩CB。全 Poisson Jacobi scheme の nilradical
についての主張ではない。本稿は新しい係数計算を行わず、受理済み actual tail、
全正次数 torsion 消滅、normalizer 比較、および新しい critical 接触模型を
組み合わせる。原点の深さと nil 三乗の零性は依然未決である。

## 1. 記号と正次数の実写像

共通 regular G19 を除いた実交叉完成環を A、その零次数 subring を A0 と書く。
正則九底は R0、共通正次数底函数は w である。既受理の相対 homogeneous 表示
は R0 上の有限 w-polynomial graded algebra A_alg を与え、その原点完成が A。
安定 tail ring を M と書くと、正確な Laurent 式は

    A_alg[w^-1] = M[w,w^-1].

これは w=0 の完成環に w=1 を代入する式ではない。先に有限 w-polynomial
algebra を使う。各 homogeneous 成分からの係数写像は A_n -> M として定まる。
受理済み T1=T2=T3=0 と degree>=4 の安定化により

    A_n -> M は同型、すべての n>=1。

全87 graph を残した同じ actual presentation の結論であり、単なる先頭項や
選択断面の推測ではない。旧49座標・全218関係との実 source binding は
M の whole nilradical が principal square-zero であることを証明している。
被約商 M_red は8次元であり、M 自体を一変数二重環や被約環へ置き換えない。

## 2. 混合積も含めた全 nil 二乗の次数

N=Nil(A) とし、N_n=N intersect A_n と書く。nilradical は連続 torus 作用で
安定な閉 ideal なので、homogeneous 成分を取っても N に入る。これは有限
次数商での weight projection と完備性から従う。有限生成 ideal は閉じており、
全成分の積の恒等式は各有限次数商を通して完成環にも移る。

a in N_i、b in N_j、i+j>=1 を取る。その M への像は両方とも nilpotent。
M の nilradical は square-zero なので像の積は零である。一方 A_(i+j)->M
は injective だから

    a*b=0 whenever i+j>=1。

特に N_positive^2 だけでなく、**N0*N_positive も零**である。したがって

    N^2=N0^2、N^3=N0^3

はすべて degree0 にある。さらに任意の正次数 c に対し、a*b*c は
a*(b*c) と書け、b*c は正次数 nil だから零となる。従って全 positive ideal
は N^2 および N^3 を annihilate する。完成に含まれる無限級数でも各次数の
係数は有限和なので、同じ主張が成り立つ。

この議論は normalized tensor 比較 kernel の全正次数が零だとは仮定しない。
使う injectivity は actual A_n->M のものである。一般点で w に殺される
kernel を原点まで一様に拡張する操作も用いない。

## 3. 三乗は比較 kernel の有限長 torsion に入る

同じ零次数 ambient 上の二つの regular normalizer を tensor した実環を E4、
K=ker(A0->E4) とする。normalization は integral で Spec 上 surjective だから、
二段の base change 後も比較写像は integral かつ Spec 上 surjective。
従ってその kernel は nil ideal である。写像は nilpotent を nilpotent に送り、
受理済み whole E4 の Nil(E4)^3=0 から

    N0^3 subset K

を得る。

conductor の逆乗法 tensor 議論により I_A0*I_B0 は K を annihilate する。
この積の支持は R0 の四つの conductor curve で尽くされ、ell は全枝上で
Noether parameter である。したがって K は k[[ell]] 上有限である。
四つの一般点での実交叉は、独立監査された以下の source 模型で与えられる。

| 一般支持 | 実交叉の nil に必要な結論 |
| --- | --- |
| CA-critical V | 新模型の nil=(x)、x^2=0 |
| CA-selfpair Pa | nil=(f)、f^2=0 |
| CB-critical V | 実 CA の外なので交叉は空 |
| CB-selfpair Qb | nil=(f)、f^3=0、f^2!=0 |

したがって N0^3 は全ての一般支持点で零となる。一変数上有限な K の
generic fiber で零となる submodule なので、N0^3 は k[[ell]]-torsion。
すなわち R0 原点に有限長で支持される。前節と合わせて

    Nil(A)^3 = Nil(A0)^3 subset tors_(k[[ell]]) K = H0_m(A0),
    length_k Nil(A)^3 < infinity,
    (positive ideal)*Nil(A)^3=0。

未知の K 全体や T0 全体を有限長としたのではない。K は一般 rank4を持つ。
有限長なのはここで特定した nil 三乗 submodule である。その長さや
maximal-ideal による annihilation 指数は計算していない。

## 4. 原点の深さおよび nil 指数との厳密な関係

既受理の exact specialized comparison を用いて L=length_k T とすれば、
先の corollary は

    depth_cone(D at P0)=20 iff L=8,
    torsion minimal generators of K = L-8

を与える。従って

    L=8  implies Nil(A)^3=0。

一方 Qb の実 family は P0 を通り、その一般点では nil の二乗が非零である。
nil ideal の二乗は coherent で支持が閉じているため、その閉包にある P0 の
局所環でも nil 二乗は非零。形式的完成と共通群 product によって失われない。
従って P0 での nil 指数は少なくとも3であり、L=8なら正確に3である。

逆向きに Nil(A)^3 が非零なら H0_m(A0) が非零なので L>8、元の cone 深さは19。
しかし Nil(A)^3=0 から L=8 は従わない。nil 三乗に属さない閉点 torsion が
残る可能性を排除していない。

共通 regular G19 product を戻すと、原点の nil 三乗は

    k[[G19]] completed_tensor_k Nil(A)^3

である。非零なら k[[G19]] 上有限自由で、その rank は未計算の transverse
有限長に等しい。位相的支持はこの B18 formal germ の19次元群方向に限られる。
projective chart では scalar を除いて18方向となる。これは局所 formal 支持の
結論であり、新しい全域 associated closure や分類・成分数は主張しない。

必要なのはこの有限長 nil-cube またはより大きな torsion K の実判定である。
本稿はその零性、元 D の全 nil 指数、原点の深さを決めず、tail を被約と
仮定せず、generic kernel の次数輸送を global annihilation に取り違えない。
実依存・独立監査の exact path と SHA256 は同名 report に束ねる。
