# 二つの楕円型 graph 成分の交叉：旗円錐と原点の冪零部分

2026-11-01。対象は被約な二成分 T39・T40 の、S61 内でのスキーム論的交叉である。
全 S61 の局所環や他成分との交叉を決定したという主張ではない。

## 1. 結論

二成分は既に証明した多項式 graph `q=Q51(u)`、`q=Q52(u)` である。
R=QQ[u₁,…,u₁₆] を重み1座標環、m=(u₁,…,u₁₆) を原点の極大イデアルとし、

\[
 I=(Q_{51}-Q_{52})\subset R
\]

と置く。添字 u₁,…,u₁₆ はこの節では16変数の一般的な略記であり、
保存済み S61 の u₁,…,u₆₁ の一部を順に選ぶ実際の対応は証明書に記録した。

primitive 一次二形式の係数写像を `A_beta:V→Λ²₀V*` とし、その5×4行列の
全2×2小行列式のイデアルを J とする。このとき

\[
 \sqrt I=J,\qquad J\text{ は素イデアル},\qquad
 J/I\simeq(R/m)(-2)^{\oplus10}.                              \tag{1}
\]

特に被約交叉は階数1一次部分の旗円錐で、次元5、次数24である。
交叉全体は原点以外ではこの被約円錐と一致する。原点に支持される冪零イデアルは
10次元、次数2に集中し、その二乗は零である。
これは直和環への分裂や「別の十点」を意味しない。

Hilbert 級数は

\[
 H_{R/J}(t)=\frac{1+11t+11t^2+t^3}{(1-t)^5},\qquad
 H_{R/I}(t)=\frac{1+11t+11t^2+t^3}{(1-t)^5}+10t^2.             \tag{2}
\]

## 2. 有限な線形代数による包含

計算は全て QQ 上で行う。I の二次部分の次元は45、J の二次部分の次元は55である。
直接の行空間計算により

\[
 I_2\subset J_2,\qquad \dim I_3=560,\qquad R_1J_2\subset I_3
\]

を確認した。両イデアルは二次式で生成されるので、これは

\[
 I\subset J,\qquad mJ\subset I                            \tag{3}
\]

を証明する。(3) から J/I は次数2だけに集中し、その次元は55−45=10。
また J⊂m² なので J²⊂mJ⊂I。この時点で J/I は平方零のイデアルである。

J が被約であることを次節で確認すると、(1) が直ちに従う。
従ってこの結果のために I の primary decomposition を探索する必要はない。

## 3. J が旗円錐の素イデアルであること

非零な階数1形式は β=ℓ⊗ω と書ける。閉性は ℓ∧ω=0 を意味するため、ω=ℓ∧m と書ける。
primitive 条件は span(ℓ,m) が symplectic 双対空間の isotropic 二平面であることを意味する。
従って階数1形式の射影化は

\[
 F(1,2;4)=\operatorname{Sp}_4/B
\]

という旗多様体である。写像は `(Cℓ⊂span(ℓ,m))↦[ℓ⊗(ℓ∧m)]`。
これは W=V(1,1) の最高重み軌道であり、対応する埋込みの最高重みは ω₁+ω₂。
その円錐 C の素イデアルを P と書くと、小行列式は C 上で零なので J⊂P。

最高重み軌道円錐の次数 n の座標空間は Cartan power `V(n,n)*` である。
この標準事実は、最高重みベクトルの n 乗が生成する既約部分表現を考えても得られる。
一次文献の記述としては Wayne A. Johnson,
[*Multivariate Hilbert series of lattice cones and homogeneous varieties*](https://bpb-us-w2.wpmucdn.com/sites.uwm.edu/dist/d/129/files/2016/04/wayne-a-johnson-1wx801b.pdf),
2015, §5.3.2, Theorem 5.3.7（印刷頁66）の一つの重みの場合を参照。
本文と直前の Cartan 積の議論を確認した。

Sp4 の Weyl 次元公式は

\[
 \dim V(a,b)=\frac{(a+1)(b+1)(a+b+2)(a+2b+3)}6
\]

なので、dim V(n,n)=(n+1)⁴。従って C の座標環の Hilbert 級数は (2) の最初の式である。

一方、J の55個の二次 Gröbner 基底を次のように完全検証した。

1. 候補の二次行空間が、小行列式の二次行空間に等しいことを QQ 線形代数で確認。
2. 先頭単項式が互いに素でない全 S-polynomial を割り算し、余り零を確認。
   互いに素な組は Buchberger の積判定を適用。
3. 得られた単項式イデアルの Hilbert 級数を、`I:x` と `(I,x)` の完全列を使う
   再帰で計算。有限次数だけの一致ではなく有理関数全体を算出。

その結果 H_{R/J} は (2) と厳密に一致する。
よって全射 R/J→R/P は全ての次数で同じ次元を持ち、核は零。
従って J=P は素イデアルであり、次元や次数だけから被約性を推測していない。

## 4. スキーム論的な読み方

二つの graph の定義イデアルは `(q−Q51(u))` と `(q−Q52(u))`。
その和による商は R/I なので、上の計算は正確に二成分のスキーム論的交叉を記述する。

完全列

\[
 0\longrightarrow \mathbb Q(-2)^{10}
 \longrightarrow R/I\longrightarrow R/J\longrightarrow0
\]

を得る。原点は embedded associated point であり、原点以外では I と J の局所化が等しい。
複素数体への基礎拡大でも同じ記述が成立する。
この交叉が全 S61 の完成局所環にどのように入るか、あるいは他の成分がそこを通るかは、
この完全列だけでは決まらない。

## 5. 証明書と再現

```sh
.venv/bin/python code/elliptic_section_intersection_20261101.py \
  --output reports/new-run/elliptic-section-intersection.json
```

主証明書は `reports/ELLIPTIC_SECTION_INTERSECTION_20261101.json`。
入力候補 `reports/ELLIPTIC_SECTION_INTERSECTION_GB_INPUT_20261101.json` は
Macaulay2 1.26.05 で探索したものだが、上の Python 検証は候補の正しさを仮定しない。
同じイデアルを生成すること、全 Buchberger 条件、Hilbert 級数を全て再確認する。
再現には Macaulay2 を必要としない。

16変数の一般 Gröbner 探索を SymPy で試すと時間を要したため、その探索を続ける代わりに、
短い候補証明書を完全検証する方式へ切り替えた。完全検証は約6秒で完了した。
