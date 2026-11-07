# smooth cubic scroll の残る A=(2,5) 二型：scalar connection への全還元

2026-11-07。主義務 G。純幾何 producer proof、別作者 review 前。
`CUBIC_SCROLL_NORMAL_ROOT_DIRECTRIX_REDUCTION.md` と
`CUBIC_SCROLL_CANONICAL_ROOT_EXCLUSION.md` の続きであり、既存 bytes は変更しない。
全 G の非優越はまだ結論しない。

## 1. 正確な入口

`S=F1⊂P4`, `H=C0+2f`, `A=2C0+5f` とする。
実既約被約 quintic `D` の rank 四 Poisson source で、`S⊂Sing D`、restriction 一般 rank 二を仮定する。
前二 note の全 primitive root の還元後に残るのは

* normal quadratic rank 二：`G=gamma ell1 ell2`, `ell_i:N→A`, `V=div gamma=C0`,
  `Delta=div(ell1 wedge ell2)=C0+3f`, `R=V+Delta`。
* normal quadratic rank 一：`G=gamma ell²`, `ell:N→A`, `div gamma=C0`。

である。rank 一では最初から第二の root 又は `R=C0+Delta` を仮定しない。
root maps は divisorial zero を持たないが、codimension 二の basepoints は許す。
`R` は reducible/nonreduced でよい。

**還元命題。** どちらの場合にも、同じ line bundle `A` へ向かう二つの独立な global rows
`ell1,ell2` を取れて、`U=S\R` の flat conormal connection に対し同一の scalar one-form `theta` を持つ。
すなわち local frames で

`dM+M beta=theta M`, `M=(ell1;ell2)`, `d theta=0`。        (1)

更に `div(det M)=R−C0`。
rank 二では二行は元の二つの normal roots、rank 一では第二行は residue の有理 primitive から構成する。

## 2. normal connection の全 multiplicity pole budget

`sigma=pi|S`、`R=div sigma` とする。Hamiltonian action を anchor `sigma#` の逆で割って
normal connection に直すと、その local matrix `beta` は

`sigma beta` が正則な rational one-form matrix

となる。従って `R` の prime `B` の multiplicity が `n` なら、`beta` の pole order は高々 `n`。
`R` の外には pole がない。この必要上界は `R` が reduced との仮定を使わない。

primitive row `ell` が parallel line を張れば
`d ell+ell beta=theta_ell ell`。
prime の generic point で row の一成分は unit なので、その成分を用いて `theta_ell` を解くと
`theta_ell` の pole も高々 `n`。
flat 性から `d theta_ell=0`。
同じ line bundle `A` の二行に対する difference は frame-independent な global rational one-form である。

## 3. rank 二：全 collision における theta1−theta2 の pole 消滅

`eta=theta1−theta2` とする。
local regular matrix `M` に対し

`beta=theta2 Id+M^(-1) diag(eta,0) M−M^(-1)dM`。          (2)

prime `B` における `m=ord_B(det M)` を取る。
`m>0` なら `M` は residue field 上 rank 一であり、二行とも非零である。
このとき

`adj(M) diag(1,0) M`

の少なくとも一成分は unit となる。
実際 `M=(a b;c d)` とするとこの積は `(da db;−ca −cb)` であり、
二行の primitivity と rank 一から全成分零にはならない。

`eta` が order `k>0` の pole を持てば、(2) の middle term は少なくとも一成分で pole order `m+k`。
一方 `M^(-1)dM` の pole order は高々 `m`、`theta2` と `beta` の pole order は §2 により高々 `n`。
従って

`k≤n−m`。                                               (3)

`m=0` でも (2) 又は §2 から同じ上界 `k≤n` がある。
全 divisor equality `R=C0+Delta` は `B≠C0` で `n=m`、`B=C0` で `n=m+1` を与える。
よって `eta` の pole は `C0` のみにあり、その order は高々一。
repeated collision の `m>1` も (3) に含まれ、generic simple collision の外を捨てていない。

closed rational one-form の simple pole は logarithmic で、その residue は定数。
一般 ample curves の residue 定理又は Gysin sequence は
`res_(C0)(eta) [C0]=0` を与える。
`[C0]≠0` なので residue は零、従って pole 自体が消える。
`H0(F1,Omega1)=0` より `eta=0`。これで rank 二の (1) が従う。

## 4. rank 一：closed residue の exact primitive

`Y=Bl_S P4`, `E=P(N)` とする。intrinsic conormal 可換性により Poisson lift は正則で、
`E` は反標準因子の multiplicity 一成分である。
inverse closed two-form の residue `alpha=res_E(omega)` は closed rational one-form。
`D'|E=p*(C0)+2Gamma` なので、pole order は horizontal prime `Gamma` で高々二、
vertical prime `p^*C0` で高々一、それ以外にはない。

generic `P1` fiber 上で `alpha` はただ一つの pole（高々二位）を持つ。
fiber residue 定理からその residue は零。従って `Gamma` 上の residue は零である。
vertical residue を `mu` とすれば、E 上の Gysin class relation は `mu p*[C0]=0`。
従って `mu=0` で、closed 性と単極性から `C0` 上の pole は消える。
従って `alpha` は全 prime で residue 零、pole は `Gamma` の二位だけである。

**有理 primitive 補題。** 有理 variety 上の closed rational one-form で、全 prime divisor の residue が零なら、
それは rational function の微分である。
ここでは `E` は rational threefold なので、function field は `C(t1,t2,t3)`。
一変数 `t3` の rational integration を `C(t1,t2)` 上で行う。
partial fraction の単極係数は対応する prime の residue で零なので、rational primitive が存在する
（必要なら有限 splitting field の primitive を trace して base field へ戻す）。
その微分を引くと `dt3` 成分が零、closed 性により他成分は `t3` に依存しない。
残る二変数・一変数へ同じ操作を帰納して `alpha=df` を得る。

`df` に `Gamma` の高々二位 pole しかないため、`f` の pole は `Gamma` の高々一位だけである。
`Gamma` の defining section を `ell` と書くと

`f=ell'/ell`, `ell'∈H0(S,N*⊗A)`。                       (4)

これは `H0(E,O_E(Gamma))=H0(S,N*⊗A)` から従う。
generic fiber 上の `alpha` は非零である。実際、`pi|S` の anchor は rank 二であり、
`pi|E` の rank 二 distribution は S の二接方向へ同型に射影し、その annihilator である `alpha` は
vertical tangent に非零でなければならない。
従って (4) は fiber 上非定数、`ell'` と `ell` は独立である。

## 5. rank 一の第二行も同じ character を持ち、determinant は R−C0

`alpha` は `pi|E` の horizontal distribution を消す。
従って `df=0` は、この distribution に沿う `f` の derivative が零、すなわち projective normal connection の
parallel function を意味する。
元の root `ell` は `d ell+ell beta=theta ell` を満たす。
(4) の ratio の水平微分を取ると

`d ell'+ell' beta=theta ell'`。

これで rank 一にも (1) が成立する。

determinant の divisor は、fiber 上の residue の実係数から求められる。
generic normal coordinates `(u,v)` を取り、blowup の fiber coordinate を `t=v/u` とする。
normal quadratic `G=gamma ell²` に対し、inverse の `du/u wedge dt` 係数、すなわち
`alpha` の vertical 成分は、固定非零定数を除いて

`sigma/(gamma ell²)`

である。これは tangent rank 二 block の Pfaffian inverse の該当成分を取るだけで、
mixed normal/tangent first jets を零と置かない。
一方 `d_fiber(ell'/ell)` の numerator は `ell wedge ell'`。
従って section の恒等式として

`det M = const * sigma/gamma`

となり、`div(det M)=R−C0`。rank 一の第二行を仮定したのではなく、closed residue から構成している。

## 6. この還元の限界

`N*(2H)=beta*Q` と `A−2H=f=L−C0` から

`H0(N*⊗A)=H0(P2,Q(1)⊗I_p)=H0(P2,T_P2⊗I_p)`。

従って残る二行は、blowdown point `p` を固定する projective linear vector fields の六次元 Lie algebra 内の
二つの独立元である。
しかし同一 scalar character を得たことから、この二 vector fields が可換又は Lie subalgebra を張るとは
まだ証明していない。その unsupported な矢印で全 source を閉じない。
rank 一・二を共通の scalar two-row frame へ還元した点だけが本 note の結論である。
全 actual quadratic bivector への extension、残る G 次元上界、cubic cone は未解決である。
