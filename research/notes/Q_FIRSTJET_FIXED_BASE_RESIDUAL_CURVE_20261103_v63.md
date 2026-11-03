# 固定 spectral base の残余曲線：二つの (2,2) 式と m 次数上界8

2026-11-03。前 packet の二 adjoint 正規形から、特性多項式を固定した後の**全三式**を導き、相対符号を除いた開 chart を**二つの明示的な Laurent 式**へ書き換えた。さらに同じ既存点の universal W₁₆ 上の一回の6×16微分により、新しい有理不変量 m が spectral base から独立と証明した。

有限な ordered-eigenvalue 拡大の後、選ばれた実像成分の残余曲線は m-line へ一般次数高々8で写る。全8点が同じ成分に属すること、正確な次数、曲線の種数、有理性・楕円性は証明していない。NEW25 の方程式でもない。

## 有限な base 座標と三つの残余式

前稿の記号を使い、Δ=x²−y²、S=x²+y²、

\[
 D=\operatorname{diag}(x,-x,y,-y),\quad
 C=\begin{pmatrix}a&1&g&h\\c&-a&i&j\\-j&h&d&1\\i&-g&f&-d\end{pmatrix}.
\]

必要な既知単元は xyΔ である。m=gj−hi、ℓ=(a+x)²−y² と置き、

\[
 d=-(ax+S)/y,\qquad c=m-\ell,\qquad f=m-x^2\ell/y^2.
\]

canonical pencil の係数を det(tC+zD)=Σpₖt⁴⁻ᵏzᵏ とする。全係数計算は

\[
 F_4=20S,\quad p_4=x^2y^2,\quad
 p_3=-2ax\Delta-2x^2S
\]

を与える。従って full characteristic base の一般点で x,y は有限な固有値選択により決まり、

\[
 a=-\frac{p_3+2x^2S}{2x\Delta}
\]

が有理的に決まる。x²,y² は Z²−SZ+p₄ の二根であり、x,y の平方根と順序の選択は高々8重の有限拡大である。ここから固定 base の係数体を

\[
 K=\mathbb Q(x,y,a,k_1,k_2),\qquad k_1=p_1,\ k_2=p_2
\]

と書く。これは full base を新たな独立パラメータで置き換えたものではなく、上記の明示的な有限 chart である。

残る四未知数は g,h,i,j。次を置く。

\[
 B_2=-4a^2x^2-2ax^3-6axy^2+3x^4+y^4.
\]

全三式は p₁=k₁、p₂=k₂、det C=Δ² で、p₂ は特に単純である。

\[
 p_2=B_2+(x+y)^2gj-(x-y)^2hi.
\]

p₁ と determinant は次のコンパクトな式で全係数を与える。

\[
 \frac{p_1}{2}
 =ad(ay+dx)+(ay+dx)(gj+hi)
   +cy(d-gh)+fx(a+hj)-xgi+yij,
\]

\[
 \det C-\Delta^2
 =(a^2+c)(d^2+f)+(gj-hi)^2+2ad(gj+hi)-\Delta^2
\]
\[
 \hspace{10mm}+cg^2+fj^2-cfh^2-i^2
       -2cdgh+2afhj-2agi+2dij.
\]

コードは直接4×4 pencil determinant とこれらの式を比較する。y²で分母を払った p₁−k₁ と det C−Δ²、および p₂−k₂ の全係数を証明書に保存した。最後は既存77項多項式と完全一致する。

## m を残した相対符号 quotient

κ=k₂−B₂ とし、

\[
 u=gj=\frac{\kappa-(x-y)^2m}{4xy},\qquad
 v=hi=\frac{\kappa-(x+y)^2m}{4xy}.
\]

すると u−v=m、(x+y)²u−(x−y)²v=κ が全係数で成り立つ。

ghuv≠0 上で、前稿の effective な相対符号 (g,h,i,j)↦−(g,h,i,j) の quotient を

\[
 R=g/h,\qquad T=gh
\]

で表す。逆には

\[
 g^2=RT,\quad h=T/g,\quad i=vg/T,\quad j=u/g.
\]

g の二つの平方根の違いはちょうど相対符号である。従って R,T は、この相対符号を除いた chart の厳密な座標であり、平方根を無断で有理関数と扱ってはいない。

次の m に依存する係数を定義する。

\[
 A_0=ad(ay+dx)+(ay+dx)(u+v)+cdy+fax,
 \qquad K_1=k_1/2-A_0,
\]
\[
 B_0=(a^2+c)(d^2+f)+m^2+2ad(u+v)-\Delta^2.
\]

三式は、この quotient 上で次の二式に正確に書き換わる。

\[
 \boxed{L_1=-cyT-xvR+fxu/R+yuv/T-K_1=0,}
\]

\[
 \boxed{\begin{aligned}
 L_2={}&B_0+cRT+fu^2/(RT)-cfT/R-v^2R/T\\
       &-2cdT+2afu/R-2avR+2duv/T=0.
 \end{aligned}}
\]

実装は j=u/g,i=v/h を代入し、全ての even Laurent monomial を

\[
 g^r h^s\longmapsto R^{(r-s)/2}T^{(r+s)/2}
\]

へ送る商写像で、p₁−k₁=2L₁、det C−Δ²=L₂ を全係数で確認した。標本による関係式の推測ではない。

## 二つの (2,2) 式、境界、交点長さ8

Q₁=RTL₁、Q₂=RTL₂ とすると、R,T についてどちらも二重次数 (2,2) の多項式である。

\[
 Q_1=-cyRT^2-xvR^2T+fxuT+yuvR-K_1RT,
\]

\[
 \begin{aligned}
 Q_2={}&cR^2T^2-cfT^2-v^2R^2+fu^2+B_0RT\\
       &-2cdRT^2+2afuT-2avR^2T+2duvR.
 \end{aligned}
\]

generic m と generic base parameters の係数体で考える。Q₂ は k₁ に依存せず、Q₁ は k₁ に線形で、その係数は −RT/2 である。Laurent ring では RT は単元である。Gauss の補題を用いると、Q₂ の非定数因子が Q₁ を割れば、k₁ の係数も割ることになり矛盾する。従ってこの generic field 上で二式は共通曲線成分を持たない。

P¹_R×P¹_T に (2,2) として同次化する。Q₁ の四つの境界辺への制限はそれぞれ非零定数倍の

\[
 fxu\,T_0T_1,\quad -xv\,T_0T_1,\quad
 yuv\,R_0R_1,\quad -cy\,R_0R_1.
\]

従って xycfuv≠0 では境界での Q₁ の零点は四隅だけである。四隅での Q₂ は順に

\[
 fu^2,\quad -cf,\quad -v^2,\quad c,
\]

で全て単元。よって共通零点は境界にはない。以上は全同次係数で検査した。

proper な交叉だから、その**スキームとしての総長さ**は

\[
 (2,2)\cdot(2,2)=8.
\]

八つの相異なる点、被約性、全八点が実際の同じ像成分に属することはここからは従わない。

## m は内在的な有理不変量で、一般 fiber で非定数

全対に対し、D が regular semisimple の開集合で

\[
 S_D=\tfrac12\operatorname{tr}D^2,\quad P_D=\det D,
\]
\[
 \boxed{m(C,D)=
 \frac{\operatorname{tr}(CD^2CD^2)-S_D\operatorname{tr}(C^2D^2)
              +P_D\operatorname{tr}C^2}
       {2(S_D^2-4P_D)}.}
\]

これは同時共役不変である。D² の二つの projector P_X,P_Y を使うと

\[
 m=-\tfrac12\operatorname{tr}(P_XCP_YC).
\]

対角形では右辺は −tr(B·J₂BᵀJ₂)/2=det B=gj−hi である。P_X↔P_Y の交換でも trace は不変なので、式に固有値の平方根は残らない。コードは full diagonal-D chart で、c,f,d の graph に依存せずこの恒等式を検査した。分子は β について12次、分母は8次で、m は4次の有理不変量である。

同じ既存 QQ 点[0] において、一回の **universal W₁₆ 上**の微分を計算した。五つの Inv₈ の微分に dm を追加した6×16行列は rank6。最初の6個の β 列の minor は QQ 上非零で、分子544 bit、分母374 bit（補助的な還元値1430）。全6微分が元の全10 Sp₄ 軌道方向を殺すことも確認した。

これは NEW25 source の rank6ではない。実際、同じ行列を保存済み actual NEW25 source 微分に合成すると rank は5のままで、dm 行は既存の五行の線形結合となる。その全係数も保存した。この一点の微分関係から、NEW25 の大域方程式を推測してはいない。

universal rank6により、m は full characteristic base に関して超越的で、実際の16次元像成分の一般 spectral fiber 上で非定数である。ordered-eigenvalue の有限拡大後もこの性質は変わらない。c,f,u,v は m の非零一次式なので、その実成分上で恒等零にはならず、上の Laurent chart はその一般 fiber に本当に交わる。

よって、その選ばれた実像成分の残余曲線は m-line へ一般有限で写り、次数は総交点長さから **高々8**。他の像成分や交点 multiplicity が8の一部を使う可能性があるので、正確に8とは言わない。この上界は明示した有限 ordered-eigenvalue 拡大と相対符号 quotient の後のものであり、全基礎体上の被覆次数へ無断で読み替えない。

## 未決事項と再現

得られたものは genus の証明ではない。残余曲線の滑らかな固有モデル、境界の正規化、実成分に属する交点の分離、正確な m-degree が未決である。一般の Sp₄ spectral plane quartic やその Prym の genus を、この残余曲線の genus と同一視しない。

NEW25 はこの一般曲線方向の未知 multisection を与えるが、その次数や方程式を今回求めていない。高次数の点評価・追加 β 点・GB は使わない。

* `code/q_firstjet_fixed_base_residual_curve_20261103_v63.py`
* `reports/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63/{report.json,certificate.json.gz}`

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_firstjet_fixed_base_residual_curve_20261103_v63.py \
  --output work/fixed-base-residual-curve-replay
```

60秒・2GB上限、実行約4.65秒。全 pencil 係数、全三式、二つの Laurent/同次式、境界恒等式、m の全勾配、universal rank6 minor と actual NEW25 source への合成を保存する。
