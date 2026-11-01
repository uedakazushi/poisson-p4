# 非零トレースの実際の零一次部分は AB で尽くされる

2026年11月1日。以下は複素数体上の被約集合の主張である。

\[
 F_0=(T01\cap\{\beta=0\})_{\mathrm{red}},\qquad
 s(Q)=Q^{12}+Q^{34}
\]

とおく。このノートは

\[
 F_0\cap\{s\ne0\}\subset\overline{AB}
\]

を示す。従来の「各層の一般点が除外される」という次元論を、全点の判定へ置き換える。
零トレースの R0 は [既に処理した](TORIC_ZERO_JET_R0_20261101.md)。
零トレースの C0、L0 の全境界は、このノートの結論には含めない。

## 1. トレースの可能性

[退化五次式の完全分解条件](TORIC_ZERO_JET_FIBER_20261101.md)により、
T01 の零一次部分では

\[
 s=-\ell_1^2-\ell_2^2,\qquad p=\operatorname{Pf}(Q)=\ell_1^2\ell_2^2,
 \qquad (J\operatorname{Hess}s)^2=0.
\]

従って非零 \(s\) は rank 一の最小冪零型、または rank 二の nil22 型に限られる。
前者では \(s=x^2/2\)、後者では \(s=xz\) にシンプレクティック変換できる。
特に \(p\) は \(s\) の支持する線形座標だけに依存しなければならない。
rank 一なら \(p\in\mathbb Cx^4\)、rank 二なら \(p\in\mathbb C[x,z]_4\)。

## 2. 最小冪零型：第三成分も二つの既知平面に戻る

旧 v05 の [全被約固定トレース還元](../../legacy/v05/proofs/minimal_nilpotent_theorem.txt)
を使う。この還元は全六十係数から出発し、三つの冪零座標を被約台で除いた後、
\(a,b,c,t_1,\ldots,t_9\) の表示を与える。六つの Jacobi 生成式を

\[
\begin{aligned}
e_1&=-12cb-2at_2-2bt_6+b, &
e_2&=5ab-2at_9+bt_8,\\
e_3&=at_3+3b^2+bt_9, &
e_4&=-8ca+2at_6+a-2bt_5,\\
e_5&=5a^2-at_8+2bt_7, &
e_6&=4at_1+12c^2-3c+4bt_4
\end{aligned}
\]

とすると、その全表示について正確な多項式恒等式

\[
\begin{aligned}
\operatorname{Pf}(Q)
={}&\frac c4x^4+\frac13(az+bw)x^3
-\frac{e_6}{12}x^4+\frac{e_4}{6}x^3z+\frac{e_1}{6}x^3w\\
&-\frac{e_5}{6}x^2z^2-\frac{e_2}{3}x^2zw-\frac{e_3}{3}x^2w^2
\end{aligned}
\]

が成り立つ。従って \(p\in\mathbb Cx^4\) は \(a=b=0\) を強制する。
残る Jacobi 式は \(3c(4c-1)=0\) だから、点は M0 または M1 の固定平面に入る。
これは第三の六次元固定トレース成分の特殊点にも適用される。
[M0・M1 の全点の判定](TORIC_ZERO_JET_MINIMAL_20261101.md)により、
T01 に属するものはすべて \(\overline{AB}\) に入る。

## 3. nil22 型：完全分解は正確に AB 正規形を強制する

旧 v06 の全固定トレース表示は十一変数
\((a,b,p_3,p_4,p_5,c,v,d,A,B,u)\) を持つ。
その全表示と九つの Jacobi 生成式は
[保存証明書](../../legacy/v06/data/nilpotent22_certificate.json)にある。
この表示は線形条件の rank49 による全六十係数の還元であり、一般点だけの正規形ではない。

九式のうち

\[
E_0=6Bb-cA+5dA-2A,\qquad E_3=6aA+Bc+3Bd-2B
\]

を用いると、全 Pfaffian の二係数は

\[
[x^3y]p=A-E_0/2,\qquad [z^3w]p=-B-E_3/2.
\]

従って支持条件だけで \(A=B=0\) となる。開枝の個々の点を調べる必要はない。
残る五生成式は

\[
\begin{aligned}
f_0&=2ub+v(c-2d+1),&
f_1&=-4ab-56uv+(c-d)(c-d+1),\\
f_2&=2av-uc,&
f_3&=a(1-2c-2d)-4u^2,\\
f_4&=b(2c-6d+3)-4v^2.
\end{aligned}
\]

閉表示の Pfaffian は

\[
\begin{aligned}
p={}&\frac b4x^4-vx^3z+\frac{d-c}{4}x^2z^2+uxz^3-\frac a4z^4\\
&+\frac{f_4}{4}x^4-2f_0x^3z-\frac{f_1}{4}x^2z^2-2f_2xz^3+\frac{f_3}{4}z^4.
\end{aligned}
\]

完全分解では \(p=h^2\) で、\(xz+2h\)、\(xz-2h\) がともに一次式の平方である。
\(h=r x^2+s_0xz+s_1z^2\) と書くと、二つの二次式の行列式の差は
\(-2s_0\)。従って

\[
s_0=0,\qquad r s_1=1/16,
\qquad p=r^2x^4+\tfrac18x^2z^2+s_1^2z^4.
\]

上の Pfaffian 恒等式と比較して

\[
u=v=0,\qquad d-c=1/2,\qquad ab=-1/16.
\]

特に \(a,b\ne0\)。\(f_3=0\) は \(c+d=1/2\) を与えるので

\[
c=0,\quad d=1/2,\quad b=-1/(16a).
\]

全二重ベクトルは正確に

\[
\begin{aligned}
Q^{12}=Q^{34}&=xz/2,& Q^{13}&=0,\\
Q^{14}&=az^2+x^2/4,& Q^{23}&=-x^2/(16a)-z^2/4,\\
Q^{24}&=p_3x^2+p_4xz+p_5z^2,\qquad &a&\ne0,
\end{aligned}
\]

となる。これは任意の磁気係数を含む AB の全固定トレース正規形である。
従って nil22 型では完全分解条件だけで AB に入る。各既約層の一般点判定や
安定化群の半連続性をこの段階で使っていない。

## 再現と依存範囲

```
.venv/bin/python code/toric_zero_jet_nonzero_trace_20261101.py --output /tmp/nonzero-trace-review.json
```

[新証明書](../../reports/TORIC_ZERO_JET_NONZERO_TRACE_20261101.json)は上の全 Pfaffian 恒等式を
多項式のまま確認し、Jacobi 生成式の明示的な組合せを保存する。
旧固定トレース表示の網羅性を使う点は依存関係として明示している。
歴史的な二十二成分の名前だけから特殊点の帰属を推論していない。
零トレースの C0、L0 内の残余を除くことが、全 \(F_0=\overline{AB}\) に向けた次の課題である。
