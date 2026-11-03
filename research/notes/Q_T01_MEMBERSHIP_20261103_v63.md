# 全Q_a族は実際のトーリック成分T01に入る

2026-11-03。全族の**帰属そのもの**を明示的な退化で証明する。
必要条件の消失、\(\Delta\ne0\) 上の公式、未飽和イデアルの特殊化からの推測ではない。

## 1. 対象と結論

\(x=(x_0,\ldots,x_4)^t\) とし、上三角の非零成分が

\[
 (Q_a)_{14}=a x_0^2,\quad (Q_a)_{23}=-\frac29x_0^2,\quad
 (Q_a)_{24}=\frac43x_0x_1,\quad
 (Q_a)_{34}=\frac23x_0x_2+x_1^2
\]

のcanonical斉次Poisson tensorを考える。T01は基本的な6パラメータ対数型族の
\(\operatorname{PGL}_5\) 飽和の既知被約閉包を表す。canonical affine coneの次元は26、
射影化した既知成分の次元は25である。

**定理。全ての \(a\in\mathbb C\) について \([Q_a]\in T01\)。**

\(a\ne0\) では原点 \([e_0]\) がrank4点であり、既存の正規化によるS61座標

\[
 u_{25}=-\frac4{3a},\qquad u_{36}=\frac3a,\qquad
 u_{38}=\frac9{2a},\qquad \text{他は零}
\]

も実際のT01に属する。\(a=0\) はS61のこのrank4正規化の外だが、
大域的な非零rank2 tensorとして同じ閉包に属する。

## 2. 可換行列による実際のT01内の普遍写像

任意の5×5行列 \(M\) に対して

\[
 W_k(M)=M^k-\frac{\operatorname{tr}(M^k)}5 I\quad(1\leq k\leq4)
\]

と置く。6個の係数 \(c_{ij}\)、\(1\leq i<j\leq4\) を使い

\[
 \Psi(M,c)=\sum_{i<j}c_{ij}(W_i(M)x)\wedge(W_j(M)x)
\]

を定義する。これはMat5×A6からcanonical Poisson空間への実際の多項式写像である。
行列は互いに可換かつtrace零なので、その線形ベクトル場も可換かつdivergence零。
従ってSchouten括弧とdivergenceは恒等的に零である。

Mが相異なる5固有値を持つ稠密開では、\(I,M,\ldots,M^4\) の独立性より、
\(W_1,\ldots,W_4\) はtracefree Cartanの基底となる。\(\Psi\) の像はその対数型族で、
既知T01に入る。T01のaffine coneは閉であり、sourceは既約だから、
\(\Psi^{-1}(T01)\) はこの稠密開を含む閉集合、すなわちsource全体である。
generic Mと全6係数を動かす像の閉包はT01そのものである。

ただし、この**実際の写像の像**がT01の全境界点を覆うとは主張しない。
T01の正規化や局所枝の尽くしも、この議論からは得ない。

標準の対数型は \(P_{ij}=b_{ij}z_i z_j\)、\(1\leq i<j\leq4\)、\(P_{0j}=0\)。
\(z_i=x_0+x_i\) への座標変更後、`code/s61_model.py` の `toric_family` の6係数族になる。
旧component catalogueのT01はedge無し・period dimension6で、そのperiod basisは全階数6。
従ってここで使う対数型は別のtoric componentの仮定ではない。

## 3. Q_aの明示的な行列表示

補助変数sについて

\[
 F(a,s)=3as^2+(9a-4)s-6=0
\]

とする。\(F(a,0)=-6\)、\(F(a,-3)=6\) なので、解では \(s(s+3)\ne0\)。
\(a\ne0\) なら二次式、\(a=0\) なら \(-4s-6\) だから、全複素aに対して解が存在する。
\(a=0\) では \(s=-3/2\) を取れる。

Nを、非零なsubdiagonal linkが順に

\[
 1,\quad s,\quad (s+3)/2,\quad -9/(s+3)
\]

である下三角nilpotent行列とする。\(N^5=0\)、
\((N^4)_{40}=-9s/2\ne0\) なのでregular nilpotentである。さらに

\[
 \lambda=-\frac{2a}{9s},\qquad
 \mu=-\frac4{9s^2(s+3)}
\]

と置くと

\[
 Q_a=\lambda(Nx)\wedge(N^4x)+\mu(N^2x)\wedge(N^3x)
 \quad\text{on }F(a,s)=0.
\]

この恒等式はEuler-wedge商だけでなく、**全150個の斉次係数**で検査した。
Fを課す前の差は正確に

\[
 \frac{F(a,s)}{3(s+3)}x_0x_1\,\partial_2\wedge\partial_4
 +\frac{F(a,s)}{6s}x_0x_2\,\partial_3\wedge\partial_4
\]

である。Nの冪はtrace零なので、右辺は第2節の \(\Psi\) の実際の値である。
これだけでも全aの帰属が従う。

Fを解いて

\[
 a=\frac{4s+6}{3s(s+3)}
\]

と書けばsourceは \(\mathbb Q[s,1/s,1/(s+3)]\)。
QQ(a)上では一般に二次拡大を使う構成であり、判別式は \(81a^2+16\) である。
この拡大上の正確な帰属は、QQ(a)上の定義イデアルの消失をも意味する。
有限標本から全aへ外挿したものではない。

## 4. canonical空間内の多項式退化も保存する

\[
 N_t=N+t\operatorname{diag}(0,1,2,3,4),\quad
 \tau_k=\operatorname{tr}(N_t^k)/5,
 \quad \widetilde W_k=N_t^k-\tau_kI
\]

を使い

\[
 \Pi_t=\lambda(\widetilde W_1x)\wedge(\widetilde W_4x)
       +\mu(\widetilde W_2x)\wedge(\widetilde W_3x)
\]

とする。\(\tau=(2t,6t^2,20t^3,354t^4/5)\)。
\(\Pi_t\) はtの多項式で、canonical divergence25係数とraw Jacobi350係数は
\(\mathbb Q(a,s)[t]\) 上で恒等零。t=0で上の150係数恒等式により \(Q_a\) へ至る。

t≠0では、固有値は \(0,t,2t,3t,4t\)。固有frameSは単位下三角で、

\[
 S_{ij}=\frac{(-1)^{i-j}\prod_{r=j}^{i-1}n_r}{(i-j)!t^{i-j}}
 \quad(i>j),\qquad S_{ii}=1
\]

（\(n_r\) は上の四link）。\(N_tS=S\operatorname{diag}(0,t,2t,3t,4t)\)、
\(\det S=1\) と逆行列を厳密に検査した。S自身にtの極があることは、
元の \(\Pi_t\) がt=0で正則であることを妨げない。

traceを引く前の同じ冪のwedgeとの差は正確にEuler場Eとのwedgeであり、
projective classは標準の対数型と一致する。全150係数でその恒等式も確認した。
対数型の6係数は

\[
 t^5\gamma_{ij},\qquad
 \gamma_{ij}=\lambda(ij^4-i^4j)+\mu(i^2j^3-i^3j^2)
 \quad(1\leq i<j\leq4)
\]

で、Pfaffianは \(288\lambda\mu t^{10}\)。generic sでは標準T01のrank4 chartに入る。
固定点 \([e_0]\) でのprojective valueのPfaffianは、元のframeでも
\(-2a/9\) でtに依存しない。従ってa≠0ではrank4正規化を保って退化を見られる。

これは、境界でCartan-limit stabilizerを見つけただけでwedge像への所属を
仮定した議論ではない。行列、wedge係数、実際のPoisson tensorまで同時に
明示的に持ち上げている。

## 5. 非零の固定一次ジェット上の直線

S61の重み付き作用 \((\beta,q)\mapsto(r\beta,r^2q)\) は、射影座標変更とtensorの
非零スカラー倍から来るのでT01を保つ。a≠0の上の点にr=aを使うと

\[
 u_{25}=-4/3,\qquad u_{36}=3,\qquad u_{38}=9a/2
\]

となる。閉包を取れば、固定した非零一次ジェット上の**全直線**

\[
 u_{25}=-4/3,\quad u_{36}=3,\quad u_{38}=c,\quad
 \text{他のS61座標}=0\qquad(c\in\mathbb C)
\]

がT01に含まれる。全ファイバーがこの直線だけとは主張しない。
c=0のS61点を、rank4正規化できない元のphysical \(Q_0\) と同一視しない。

## 6. 一般被覆・飽和レシピとの関係と検証範囲

[全域の飽和素イデアル表示](TORIC_AFFINE_MODIFICATION_20261101.md)は
\((f,\widehat\Delta(q-\bar q)-F_i):\widehat\Delta^\infty\) の縮約を使う。
今回その未飽和式を \(\Delta=0\) へ代入して十分条件とはしなかった。
全16変数の巨大な飽和GBも計算していない。

探索で使った既証の全域必要条件274本は、正しいclosed jet座標で全て \(Q_a\) 上で零。
多項式根座標 \(X_{\rm pol}\) も零で、canonical一次ジェットは
\(c_4=c_6=\Delta=0\) だった。これらは帰属の証明には使わない。
前節までの実際の多項式写像と退化が十分性を与える。

主計算：[q_t01_membership_20261103_v63.py](../../code/q_t01_membership_20261103_v63.py)。
別の必要条件点検：[q_t01_global_constraint_checks_20261103_v63.py](../../code/q_t01_global_constraint_checks_20261103_v63.py)。
保存先：[Q_T01_MEMBERSHIP_20261103_v63](../../reports/Q_T01_MEMBERSHIP_20261103_v63/)。

```sh
python code/q_t01_membership_20261103_v63.py --output work/q-t01-fresh
python code/q_t01_global_constraint_checks_20261103_v63.py --output work/q-t01-constraints-fresh
```

主証明書は約1秒、全274必要条件等の別点検は4.5秒。
全T01の局所正規化、枝の個数、全Q局所モデルの尽くしはここでは結論しない。
