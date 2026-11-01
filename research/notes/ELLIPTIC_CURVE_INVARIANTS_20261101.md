# 一次ジェットの canonical Pfaffian curve の厳密不変量

2026-11-01。コード `code/elliptic_curve_invariants_20261101.py`、証明書
`reports/ELLIPTIC_CURVE_INVARIANTS_20261101.json`。

LG(2,4) の primitive Plücker 行列を Z とし、
\(a_\beta(Z)_i=\sum_{j<k}\beta_{i,jk}Z_{jk}\) とする。canonical curve の行列

\[
 \Phi_\beta=\begin{pmatrix}\star Z&a_\beta(Z)\\-a_\beta(Z)^t&0\end{pmatrix}
\]

に Fisher–Sadek,
[On genus one curves of degree 5 with square-free discriminant](https://arxiv.org/html/1406.6953v1),
§2 の Ω の cyclic derivative construction と Theorem 2.3 を適用する。
実装は二次・三次多項式の微分を定数係数 tensor の縮約に変換し、整数演算で評価する。
有理 β の分母は第五行・列の congruence scaling で払ってから戻す。
c₄,c₆,Δ の β 次数はそれぞれ8,12,24。

正規化は Fisher,
[The Hessian of a genus one curve](https://arxiv.org/html/math/0610403v2),
§3–4 の Hesse 五次模型 u₅(a,b) によって独立に確認した。
\(D=ab(a^{10}-11a^5b^5-b^{10})\)、c₄ の Hessian 公式、c₆ の Jacobian 公式を
使い、五つの整数 (a,b) で三不変量すべてが一致し、Δ=D⁵ となることを確認した。
β の3倍によるスケーリングも次数と整合する。

両楕円族について λ=2,3、p=(0,0,0,0),(1,2,3,4),(2,3,4,5) の計12例を正確に計算した。
全例で Δ≠0 なので、上記一次資料の Theorem 2.1 により表示された Pfaffian curve は滑らか。
各 λ について c₄,c₆ そのものが p によらず、次の式と一致した：

\[
 (c_4,c_6)_{q(5,1)}=
 \left((4/\lambda^3)^4H_4(\lambda,-1),(4/\lambda^3)^6H_6(\lambda,-1)\right),
\]
\[
 (c_4,c_6)_{q(5,2)}=
 \left(\lambda^{-12}H_4(1,\lambda),\lambda^{-18}H_6(1,\lambda)\right).
\]

**これらの全 λ,p にわたる恒等式は、この有限標本証明書だけからは主張しない。**
general p の構造的復元は別の kernel / Grassmannian 議論で扱う。
λ=2 の全六例では j=6761990971/5153632 で、source Hesse curve の j と一致する。

低次不変量による短縮も調べたが、c₄ は span{I₄²,tr(A⁴),vᵀJAv} に入らず、c₆ も
span{I₄³,I₄tr(A⁴),I₄vᵀJAv,vᵀJA³v,T(v,v,Av)} に入らない。
12個の整数 primitive β による exact coefficient rank はそれぞれ 3→4、5→6。
これは該当する短い候補公式を排除する反例であり、全不変式環を計算した結果ではない。
探索データは ignored scratch `work/canonical_invariant_fit.json` にある。

トーリック inverse algebra の特定の巡回生成元 y の三次式を y²=f(z) と見なした曲線は、
canonical curve と同じ j を持たない。β* では前者が1728、後者が
7980244148583563362683690248/2081757871309100985 である。
ただし三つの解の étale algebra が canonical curve の非零二捩れ点の代数に同型である
可能性はこれで排除されない。三次代数の任意の生成元は、楕円曲線の二分点の
Weierstrass x 座標とは限らないためである。

```sh
.venv/bin/python code/elliptic_curve_invariants_20261101.py \
  --output /tmp/elliptic-curve-invariants-review.json
```
