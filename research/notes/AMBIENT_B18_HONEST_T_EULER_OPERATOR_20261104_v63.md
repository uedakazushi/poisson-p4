# 実 T の完全116生成・936関係：形式 Euler cycle operator

2026-11-04、v63。比較 kernel の深さ判定へ向け、actual 特殊比較 module

\[
 T=\ker\bar V/\operatorname{im}\bar X
\]

の**完全な有限形式 operator presentation**を構成した。
116生成、936関係、74,176 typed nodes、70,017 coefficient slots、証拠3.99 MBである。
係数を348環へ展開した行列でも、(T\) の長さを求めた結果でもない。

## 1. 実 source の768 cycle

\(\bar V: \bar R^{56}\to\bar R^3\) の source は、既存 (Q48\) と (-U8\) の順であり、
target は実 \((\tau_A,\tau_B,\tau_A\tau_B)\) である。
\(\bar X\) の最初640列を

\[
 \begin{pmatrix}B48\\P8\end{pmatrix}
\]

とし、後の128列を

\[
 \begin{pmatrix}\operatorname{Syz}Q48\\0_8\end{pmatrix}
\]

と定めた。除いた (P0\) は actual (UP=QB\) の定数行により正確に整合する。
128列は既受理の actual Res4 operator
\(\operatorname{Res}_4(N^t)\operatorname{Res}_4(S_{12\times32})\) の実 coefficient entry であり、
native source の符号を変更していない。

全 source column と coefficient callback を node ID・file hash で保存した。
元の正則 (r/g\) 関数はその実関数のまま保持し、特殊化した分母を割り直さない。
全関数の \(\ell=0\) 制限は、CA F9 の unit 座標を使う正確な
\((z_A,r_A,s_A)=(z_A,r,r)\) の formal composition である。

actual (UP=QB\) と (Q\operatorname{Syz}Q=0\) により、
この768列は全て \(\bar V\) の実 cycle である。

## 2. unit column change と三つの Koszul sector

先行 packet の実正則 (m:24\times56\) と選択24列 (M\) は、

\[
 d_1m=\bar V,\qquad n=E_SM^{-1},\qquad mn=I_{24}
\]

を満たす。未選択の32 source 行を (J\) とすると、actual unit column change は

\[
 \bar V\longmapsto[d_1,0_{32}],\qquad
 \bar X_j\longmapsto(m\bar X_j,\bar X_{J,j})
\]

である。従って \(\ker\bar V\) は実 source 同型で
\(\operatorname{im}d_2\oplus\bar R^{32}\) と同定される。

実 (d_1,d_2,d_3\) は三 regular sequence

\[
 (z_7,r^3),\quad(z_7+c r^2,r^3),\quad(z_7,r^2)
\]

の literal Koszul direct sum である。

## 3. corrected Euler contraction

各 sector を adapted 座標 (u\) に移す。
中央 sector は (z=u-c r^2\)、他は (z=u\) であり、(r\) は固定する。

\[
 g=(u_1,\ldots,u_7,r^n),\quad
 E=\sum_{i=1}^7u_i\partial_{u_i},\quad
 \sum u_i v_i+r^n v_8=0
\]

に対して、Koszul の符号を (d_2(e_i\wedge e_j)=g_i e_j-g_j e_i\) とし、

\[
 w_{i8}=(E+1)^{-1}\partial_i v_8,\qquad
 w_{ij}=(E+2)^{-1}(\partial_i v_j-\partial_j v_i),\quad i<j\le7
\]

と定める。**Euler inverse は微分後の出力係数に作用する。**
出力の (u\)-次数 (d\) に対して (1/(d+1)\)、(1/(d+2)\) を掛けるだけであり、
正の整数以外をこの contraction で割らない。

第一の七成分について

\[
 \sum_i u_i(\partial_i v_j-\partial_j v_i)
 =(E+1)v_j+r^n\partial_jv_8
\]

なので、Euler inverse と (u_i\) の次数差を正しく扱うと

\[
 \sum_i u_iw_{ij}=v_j+r^n w_{j8}.
\]

Koszul differential は最後の項を引く。
第八成分は

\[
 \sum_i u_iw_{i8}=v_8-v_8|_{u=0}=v_8
\]

である。最後の等号は cycle を (u=0\) へ制限した
\(r^n v_8(0,r)=0\) と domain \(\mathbb Q[[r]]\) から従う。
従って全形式級数について (d_2w=v\) である。
基本 cycle \((-u_2,u_1)\) は (w_{12}=1\)、
\((-r^n,0,\ldots,u_1)\) は (w_{18}=1\) を与える。
全28 wedge・2128個の多項式倍の固定 fixture でも符号と出力次数を確認した。
有限 fixture を全次数の証明の代用にはしていない。

この operator を actual (m\bar X_j\) に適用し、元の (z\) 座標へ戻す。
**先に (I^2\) や348環の normal form へ落として微分してはいけない。**
Euler inverse は一般に algebraic function とは限らない形式係数を作るが、
完成環 \(\bar R\) の正則形式級数として定義される。
contraction 自体は \(\bar R\)-線形ではない。
768個の各 coefficient column を選んだ後、その列から自由加群写像を \(\bar R\)-線形に延長する。

## 4. 完全な T presentation

実 generator は (G_2\) の84 wedge と未選択 source32座標の計116個。
関係は、実 (d_3\) の168列と

\[
 \bigl(w_j,\bar X_{J,j}\bigr),\qquad j=1,\ldots,768
\]

の計936列である。各 (w_j\) は上の operator により (d_2w_j=m\bar X_j\) を満たす。
\(\ker d_2=\operatorname{im}d_3\) は既受理の actual literal resolution で完全なので、
これは余分な未確認 source syzygy を落としていない、実 (T\) の完全な cokernel である。
Syz(P\) を計算したという主張ではない。

`actual-honest-T116-936-Euler-operator.json` は、全 sourceX、原 source の768 cycle 入力、
84×768 coefficient readback、literal (d_3\)、全116×936 operator columns を保存する。
大きな (mX\) の scalar 展開を行わず、既知 matrix product の係数 callback として保存する。
callback は実 matrix entry の有限和であり、新しい未定係数を導入していない。

## 5. 未計算の係数・長さ

この完全 presentation ができたので、正確な係数 restriction を用いる次の段階を定義できる。
現在、(X(0)\)、(X\) の一次、unit minimalization、348環の係数、
\(\operatorname{length}T\)、\(K\) の torsion・元の深さは計算していない。
coefficient の derivative/Euler は必ず元の実形式関数に作用させる。
既存の非単位 quotient に対する必要な source Taylor 情報は個別に audit する。
単なる特殊化 kernel や一次情報から長さを推測しない。

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_honest_T_Euler_operator_20261104_v63.py
```
