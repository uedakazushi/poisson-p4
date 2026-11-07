# cubic scroll blowup の bivector は exceptional restriction から一意

2026-11-07。主義務 G。純 bundle 証明、独立 review 前。
`CUBIC_SCROLL_QUADRIC_NET_GEOMETRY.md` の全 incidence 表示を用いる。
新 CAS、sample、minor、jet は必要ない。

## 1. 全域の injectivity

`Y=Bl_S P4=P(F) --p--> P2`、`H=O_P(F)(1)`、`B=p*O(1)` とし、
`0→F→O^5→O(1)^2→0`、`E=2H−B` とする。

**補題。** `H0(Y,Λ²T_Y(−E))=0`。
従って regular bivectors の restriction

`H0(Y,Λ²T_Y) → H0(E,(Λ²T_Y)|E)`

は injective である。特に E を Poisson divisor とする二つの bivectors は、
intrinsic `pi_E∈H0(E,Λ²T_E)` が等しければ Y 全体で等しい。

**証明。** `p` は全 base 上の smooth P2-bundle なので relative tangent sequence

`0→T_(Y/P2)→T_Y→p*T_P2→0`

は vector bundles の完全列である。
その exterior-square filtration を `−E=−2H+B` で twist した三つの graded pieces は

1. `Λ²T_(Y/P2)(−E)=O_Y(H−B)`、
2. `T_(Y/P2)(−2H+B)⊗p*T_P2`、
3. `p*Λ²T_P2⊗O_Y(−2H+B)`。

最初の等式は `K_(Y/P2)=−3H+2B` から従う。
後二つを任意の fiber P2 へ制限すると、固定 vector space との tensor を除いて
それぞれ `T_P2(−2)`、`O_P2(−2)`。
Euler sequence と `H1(P2,O(−2))=0` によりどちらも H0 は零。
従って後二 graded pieces の全 global H0 は零である。

最初の H0 は `H0(P2,F*(−1))`。
dual Steiner sequence を twist すると

`0→O(−2)^2→O(−1)^5→F*(−1)→0`

であり、全 H0 と左の H1 が零なので、これも零となる。
filtration の全 pieces の H0 が零だから補題が従う。□

この proof は formal first jet の一意性を ambient global section の一意性へ推測したものではない。
必要な vanishing を全 Y 上で直接証明した。

## 2. 残 A25 source の parameter count への正確な使い方

`CUBIC_SCROLL_REMAINING_ROOT_SCALAR_CONNECTION.md` の same-A two-row frame の還元が成立するとする。
同 note 自体は別に独立 review が必要であり、本補題から逆向きに仮定しない。
`M=(ell1;ell2)`、`A=2C0+5f`、`div detM=R−C0` とする。
`sigma=pi|S` は

`sigma=nonzero_constant * s_(C0) detM`

である。global divisor が同じ二 section の比は proper S 上の nonzero constant だからである。
E 上の rank 二 distribution は `alpha` の kernel であり、一般点で S へ同型に射影する。
従って `sigma` と `ker alpha` は `pi_E` を一意に定める。
前節の injectivity により、これらは全 Y の pi、従って元の P4 の pi も一意に定める。

node では同じ A の二 roots と `V=C0` に対し、horizontal residues は `lambda,−lambda`、
vertical residue は class relation から零。
従って `alpha−lambda dlog(ell1/ell2)` は全 divisorial poles を失い、
proper rational E 上の正則一形式として零となる。
つまり `ker alpha=ker dlog(ell1/ell2)`。
各 root の independent nonzero scale は kernel を変えず、detM と sigma を全体 scalar 倍するだけ。
従って **projective pi は unordered pair of projective global rows** だけから高々一つである。
scalar theta や ambient extension coefficients を追加 parameter として数える必要は無い。

cusp では `alpha=d(ell2/ell1)` とし、ell1 が元の normal root、ell2 は residue の primitive から得た第二行。
第二行に divisorial zero があり得ることは許す。
`ell1→a ell1`, `ell2→b ell2+c ell1`、`ab≠0` は kernel を変えず、sigma も scalar 倍するだけ。
この三次元 affine-triangular presentation fiber を、全 admissible source を保つ場合だけ引ける。
任意の二行の GL2 mixing が node の dlog ratio を保つとの誤った商は使わない。

## 3. 条件付きの最弱 dimension bridge

固定 S に対する root-pair presentation locus の dimension が、上の真の presentation fibers を引いた後で
六以下と証明されたなら、全 PGL5 飽和 projective 像は

`24+6−dim Aut(S)=24+6−6=24`

以下となる。Aut(S) は source と exceptional data を保つ全六次元群であり、
`PGL5×source` 上の対角作用は第一因子で free なので、source 点の stabilizer が増大してもこの六次元 fiber を失わない。

例えば node の二つの P2 vector fields が、点 p の tangent plane 上でそれぞれ nonzero nilpotent 二×二 matrix を持つとする。
各 nilpotent cone は dimension 二、二つの lower row は各二係数。
二つの root scales を引けば `2+2+2+2−2=6`。
ただし **全 source がこの nil/nil locus へ入ることは本 note では証明していない**。
semisimple/scalar/zero tangent matrices、cusp、全 quartic partitions と full-P2 fibers を勝手に除外しない。
この補題は、その全境界分類が完成した場合に全 global image へ移す正確な bridge だけを与える。
