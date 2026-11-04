# 元の NEW25 族から六変数 β・N₁₂・m の有理回路へ

2026-11-03。既存の三変数 Poisson 族と評価点を、既に証明済みの有限同値を使って整理する。**六つの source 変数で primitive β と全二共変量、N₁₂,m を与える1025-nodeの厳密な有理回路**を保存した。高次数 fitting、追加点、24次評価は使わない。

六変数から五変数の rational quotient が得られたとは主張しない。残る一次元の first-jet source redundancy の積分は未決であり、後述する小さな raw-jet tangent kernel が次の具体的問題である。

## 証明済みの二つの source gauge

元の三変数族は `R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63` に全次数で固定されている。n₁≠0 上で、そこで証明済みの有限変換 a=1/n₁,ε=1 により

\[
 (n_0,n_1,n_2)\longmapsto (u,1,v),\qquad
 u=n_0,\quad v=n_1n_2
\]

とできる。元の tensor scalar を c_old、評価点を (pₓ,pᵧ,p_z,p_w) とすると、同時に

\[
 c=c_{\rm old}/n_1^2,\quad
 X=p_x/n_1,\ Y=p_y/n_1^2,\ W=p_w/n_1,
\]
\[
 Z_0=p_z/n_1^3+(n_1^{-3}-1)/(6u)
\]

へ移る。scalar を省くと同じ normalized first-jet invariant の比較にならない。

次に、族を実際に保存する既証の modular Ga flow を使う。u≠0 上で時刻 −3Y/(4u) を選ぶと評価点の y 座標が0になり、

\[
 Z=Z_0-3XY/4-3vY^2/8-Y/(64u)
\]

となる。この二つは full tensor と評価点の有限 projective 同値であり、primitive 正規化後の違いは Sp₄ である。従って Sp₄ invariant を失わず、一般 source をこの chart で扱える。

## 六変数の raw value と linear jet

残る変数を u,v,c,X,W,B とする。最後の B は評価点での Pᵧ_z 成分で、

\[
 B=-4uZ/3-W/1536+W^2-X^2/64-2/9.
\]

u≠0 なのでこれは Z の代わりの有理座標である。逆式は

\[
 Z=-\frac3{4u}(B+W/1536-W^2+X^2/64+2/9).
\]

A=1/36−256W/3 とすると、scalar を掛ける前の値は

\[
 C_0=\begin{pmatrix}
 0&A&3XA/4&8/9\\
 -A&0&B&4X/3\\
 -3XA/4&-B&0&X^2\\
 -8/9&-4X/3&-X^2&0
 \end{pmatrix},\qquad \operatorname{Pf}C_0=8B/9.
\]

linear jet は、保存済みの元の全 P(u,1,v) を x,y,z,w で微分し、評価点 (X,0,Z,W) に代入して作る。Z は定数項の B へだけ現れ、linear jet 自体は u,v,X,W のみで決まる。コードはこの全値・全一次係数を元の族から直接再構成する。非退化値の条件は cB≠0である。

この chart の一般部には uv(32u+v)cB≠0を取る。最初の条件は元の index3・滑らかな Q の一般部と一致する。回路が境界に延長できても、その境界に同じ generic component の議論を自動的に適用しない。

## 分母 cB だけの rational Darboux frame

行列

\[
 R=\begin{pmatrix}
 1&0&0&0\\
 0&0&0&9/(8c)\\
 -3X/2&1&0&-9A/8\\
 -9X^2/(8cB)&0&1/(cB)&-27XA/(32cB)
 \end{pmatrix}
\]

を取る。全六変数を残して

\[
 R(cC_0)R^{\mathsf T}=J
\]

を確認した。pivot は Pₓ_w=8/9で、A の非零性を仮定しない。平方根を使わない。

raw linear jet L と新しい変数 q に対し

\[
 L_R(q)=cR\,L(R^{-1}q)R^{\mathsf T}
\]

を作り、従来通り τ=(L_R)₀₁+(L_R)₂₃、α=dτ/5、v_α=Jα により

\[
 L_{\rm prim}=L_R-2(\alpha\cdot q)J+qv_\alpha^{\mathsf T}-v_\alpha q^{\mathsf T},
 \qquad \beta=-JL_{\rm prim}J
\]

とする。全 primitive trace、全 dβ=0、元24係数への基底復元を全係数で検査し、受理済み16基底の β₁₆ を得た。計算された全有理式と frame、primitive perspective α を保存した。

## N₁₂ と m を元の族へ厳密に結び付ける

前 packet で全係数を検証した二共変量 C,D を、β₁₆ の二次式として回路に合成する。ここで C は上の raw value C₀ とは別の adjoint 共変量である。回路は次を全て出力する。

\[
 F_4=\operatorname{tr}(-C+2D)^2,\quad
 S_D=\operatorname{tr}D^2/2,\quad P_D=\det D,
\]
\[
 \delta_8=S_D^2-4P_D,
\]
\[
 \boxed{N_{12}=\operatorname{tr}(CD^2CD^2)
             -S_D\operatorname{tr}(C^2D^2)+P_D\operatorname{tr}C^2,}
\]
\[
 \boxed{m=N_{12}/(2\delta_8).}
\]

さらに五つの Inv₈ trace 基底と det(tC+zD) の全5係数を出力する。P_D は4×4 Newton 恒等式で作る。N₁₂を六変数で全面展開せず、matrix multiplication の有限 arithmetic circuit として保存する。β 自体の有理式、元の族への有限 frame と全係数復元を既に検査したため、これは標本一致による source binding ではない。

β の定義域は上記 cB の開、m には追加で δ₈≠0が必要である。分母を恒等的な単元と取り違えない。

## scalar の重みと残る quotient 問題

R(c)=diag(1,c⁻¹,1,c⁻¹)R(1) である。この式と primitive 化の equivariance により、Cartan weight (w₁,w₂) の β₁₆ 座標の c 指数は

\[
 (1-w_1-w_2)/2
\]

となる。コードは16座標全ての有理恒等式として確認した。従って invariant の c 指数は β 次数の半分となり、回路の全 node の grading でも

\[
 F_4\sim c^2,\quad\operatorname{Inv}_8,\delta_8\sim c^4,
 \quad N_{12}\sim c^6,\quad m\sim c^2
\]

を検査した。これは scalar の自由度を除いたという意味ではない。affine quotient にはこの自由度も含まれる。

一般 NEW25 quotient は5次元で、この六変数 source には一般に一次元の冗長性が残る。その first-jet に限る最後の同値を既存 full tensor の Ga と混同しない。N₁₂が24次式に二次で現れる可能性を仮定せず、まだ multisection 方程式を得てはいない。

次の bounded 問題は、raw value＋linear jet の30行に、六つの source 方向、全GL₄の16方向、perspective4方向を並べた30×26 matrix の rational kernel である。c∂c を scalar 列に選び、全行を c で割ると、係数は u,v,X,W,B の五変数多項式となる。従来の実 rank25証明と有限 gauge の同値性は期待する rank を説明するが、具体的な kernel と非零 source 部分は別に証明する。kernel vector field の積分や五つの有理 first integrals は、その後の課題である。

## 保存と再現

* `code/q_new_family_six_parameter_covariants_20261103_v63.py`
* `reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/{report.json,certificate.json,source-six-circuit.json.gz}`

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_six_parameter_covariants_20261103_v63.py \
  --output work/new25-source-six-replay
```

60秒・2GB上限、実行約0.70秒、1025 nodes。点の評価は一度も行わず、全 raw jet・frame・primitive identity と scalar grading を記号的に検査する。新しい五変数 quotient、integrability、NEW25 の新しい方程式は主張しない。
