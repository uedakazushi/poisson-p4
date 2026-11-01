# 五次元可換代数模型の境界：20の代数型とPoisson写像

2026年11月1日。半単純部分だけに制限していたv28の模型を、非半単純代数にも適用する。
ここで20と数えるのは**五次元可換単位的代数の同型類**であり、Poisson成分やPGL5軌道ではない。

## 1. 全代数模型の意味

V=C⁴、Jを固定した非退化反対称行列とする。
可換な積 m:Sym²V→V に tr(m_u)=0 を課し、
q(u,v)=tr(m_u m_v)/3 と置く。v28で証明した96二次式は

\[
A=\mathbb C1\oplus V,\qquad
(a,u)(b,v)=(ab+q(u,v),av+bu+m(u,v))
\]

の結合性にスキーム論的に同値である。
逆に任意の五次元単位的可換代数で V=ker Tr とするとこの形になる。
実際、Tr1=5なので分解は標準的で、Tr(L_uL_v)=5q(u,v) とブロックの直接計算から
3q(u,v)=tr(m_u m_v) が従う。したがってβを固定しない36変数の全模型Xは、
単位とtrace-zero基底を固定した五次元代数のスキームそのものである。

五次元可換代数の基底付きスキームは既知の結果により複素数体上で既約である
（Mazzola, *Generic finite schemes and Hochschild cocycles*, 1980, Corollary4。
Poonen, [*The moduli space of commutative algebras of finite rank*](https://math.mit.edu/~poonen/papers/moduli.pdf),
Remark1.2, §6にも適用範囲が明記される）。単位を最初の基底に固定する操作は連結群の束であり、
さらに残る四基底をtrace-zeroにする操作は、各ベクトルからTr/5倍の単位を引くものなので、
単位固定模型は A⁴×X と同型になる。
従ってXの被約台は既約16次元で、半単純開部分 GL4/S5 が稠密である。
ここではX全体の被約性・正規性を別途主張しない。

## 2. 非半単純部分にも定義されるPoisson写像

A上の正則表現で、Vの基底の乗法行列を

\[
L_i=\begin{pmatrix}0&q_i\\ e_i&m_i\end{pmatrix}
\]

とする。結合性によりこれらは可換なので、P(A)上の対応する射影ベクトル場も可換である。
その四つのベクトル場をJで組み合わせた双ベクトルはPoissonである。
1+yのアフィンチャートでは、その係数行列は

\[
F(y)=I+m_y-yq_y,\qquad P(y)=F(y)JF(y)^t.
\]

これは非半単純代数でも成立し、P(0)=Jを満たす。
四次項は反対称性で消え、残る三次項はradialなので、射影四空間の二次双ベクトルを与える。
元の61座標への正規化は既存adapterで行う。

この写像Φ:X→S61は多項式写像で、稠密な半単純部分では既存のtoric写像と一致する。
従って**非半単純部分の像もT01に含まれる**。
ただしΦの像がT01の全閉点を覆うことは、稠密性だけからは従わない。
以下は代数模型の境界を網羅するものであり、T01境界全体の網羅性とは区別する。

PfaffianはSchur補行列から

\[
\operatorname{Pf}P(y)=\det F(y)=N_A(1+y).
\]

局所積分解 A=∏A_α、dim A_α=n_α、剰余characterをχ_αとすると、

\[
N_A(1+y)=\prod_\alpha(1+\chi_\alpha(y))^{n_\alpha}.
\]

従って実際の退化五次因子の超平面の重複度は局所代数の長さで与えられる。
全体がlocalならtrace-zero元はnilpotentなのでアフィンノルムは1、射影因子は5Hとなる。

## 3. 代数型の尽くしと検証

用いた一次資料はPoonen,
[*Isomorphism types of commutative algebras of finite rank over an algebraically closed field*](https://math.mit.edu/~poonen/papers/dimension6.pdf),
Table1および§1（2008年論文の2025年訂正版）である。
標数零のlocal型は長さ1,…,5について1,1,2,4,9個。
長さの和が5となる直積を取ると20型になる。
元の代数分類は文献の結果であり、この研究の新規な分類として扱わない。

コードはTable1の標数零の商イデアルを入力し、標準単項式、正則乗法、結合性を計算する。
各型で一つの明示的なtrace-zero基底変更を行い、元のS61の全134式、150係数の復元、
独立な射影Jacobi恒等式を検証した。ノルムの因数分解、61座標、正則乗法行列も保存した。

Aut(A)の次元は全Leibniz条件によるDer(A)の次元として計算した。
標数零なので代数群Aut(A)は滑らかであり、GL4代数軌道の次元は16−dim Aut(A)となる。
下表の写像微分階数は、選んだ点での厳密な有理行列階数である。
軌道次元と一致した行では像の次元も確定する。一致しない行では、その値だけで
一般階数や像次元を決定しない。

| 代数型 | 代数軌道次元 | 検証点での写像階数 |
|---|---:|---:|
| A1 x A1 x A1 x A1 x A1 | 16 | 16 |
| A1 x A1 x A1 x A2 | 15 | 15 |
| A1 x A1 x A3_curv | 14 | 14 |
| A1 x A1 x A3_sq | 12 | 12 |
| A1 x A2 x A2 | 14 | 14 |
| A1 x A4_curv | 13 | 13 |
| A1 x A4_21_r1 | 11 | 11 |
| A1 x A4_21_r2 | 12 | 12 |
| A1 x A4_sq | 7 | 6 |
| A2 x A3_curv | 13 | 13 |
| A2 x A3_sq | 11 | 11 |
| A5_curv | 12 | 12 |
| A5_211_s2 | 10 | 10 |
| A5_211_s1 | 11 | 11 |
| A5_22_a | 10 | 10 |
| A5_22_b | 9 | 9 |
| A5_31_r1 | 6 | 5 |
| A5_31_r2 | 8 | 8 |
| A5_31_r3 | 9 | 8 |
| A5_sq | 0 | 0 |

型名の正確な商イデアルは `code/toric_algebra_boundary.py` の `local_catalog` にある。
微分はGL4の16方向で取り、Euler垂直成分を除く。
primitive一次部分の条件は全族で保たれ、perspective方向との交わりが零なので、
この階数はS61正規化後の階数と等しい。

## 再現・残余

```sh
.venv/bin/python code/toric_algebra_boundary.py --output reports/new-toric-algebra-boundary.json
```

証拠は `reports/TORIC_ALGEBRA_BOUNDARY_20261101.json`。
全20型はそれぞれ一般のGL4基底変更を許す族としてPoisson構造を与えるが、
異なる代数型のPoisson像は一致し得る。源の代数同型類をPoisson同型類として数えない。
次は三つの階数不足の族の実際のファイバーと像の同定、そしてこの代数模型から
漏れる可能性があるtoric境界点の制御である。
