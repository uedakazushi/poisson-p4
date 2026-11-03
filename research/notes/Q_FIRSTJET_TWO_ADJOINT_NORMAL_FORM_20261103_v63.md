# 二 adjoint の七座標 hypersurface 切片

2026-11-03。四つの普遍関係式を、regular semisimple D と二つの torus gauge の開集合上で厳密に消去する。**全 W の共変量対の像を含む四式の零点が、七座標・一方程式へ書き換わる。** 同じ既存 NEW25 点により、この chart が実際の像成分と NEW25 の両方に交わることも証明した。これは NEW25 の方程式ではなく、その認識問題を置く周囲の6次元切片である。

## 行列、開条件、torus gauge

J=diag(J₂,J₂) の規約で

\[
 D=\operatorname{diag}(x,-x,y,-y),\qquad
 C=\begin{pmatrix}
 a&b&g&h\\ c&-a&i&j\\-j&h&d&e\\i&-g&f&-d
 \end{pmatrix}.
\]

これらは sp₄ の全一般形である。D をこの形にするためには、その固有値を順序付きで選ぶ有限拡大が許される。以下は

\[
 xy(x^2-y^2)be\ne0
\]

上の計算である。D の regular 条件と b,e の gauge 条件を省かない。

D を保つ torus T=diag(t,t⁻¹,u,u⁻¹) による共役では b↦t²b、e↦u²e。b,e が単元だから t²=b⁻¹、u²=e⁻¹ の有限 étale 拡大で b=e=1 とできる。四つの平方根のうち同時符号反転は中心−Iで C を変えず、effective な残りの相対符号は

\[
 (g,h,i,j)\longmapsto(-g,-h,-i,-j)
\]

である。D の固有値順序を変える有限 Weyl の選択も残る。ここで有理的な gauge 選択や quotient の有理性を仮定してはいない。

## 三座標の線形消去

前 packet の q₁=tr C²,q₂=tr(CD),q₃=tr D² と四式

\[
 f_1=q_1-2q_3,\quad f_2=q_2+q_3,
\]
\[
 f_3=\operatorname{tr}C^4+4\operatorname{tr}D^4-3q_3^2,
\]
\[
 f_4=-\operatorname{tr}C^4+6\operatorname{tr}(C^2D^2)
       +12\operatorname{tr}(CD^3)+8\operatorname{tr}D^4
\]

を使う。Δ=x²−y²、m=gj−hi、L=(a+x)²−y² と置く。

f₂=0 は d を線形に決める。b=e=1 の下で f₁ と H=(f₃+f₄)/12 は c,f について線形であり、その係数行列式は y²−x² である。解は

\[
 \boxed{d_0=-\frac{ax+x^2+y^2}{y},\qquad
 c_0=m-L,\qquad f_0=m-\frac{x^2}{y^2}L.}
\]

Δ は最終的な graph の式では消えるが、消去の同値性には Δ の可逆性が必要なままである。Δ=0 へ同じ結論を延長していない。

実装は行列 trace からこれらを直接検査した。さらに z₁=d−d₀,z₂=c−c₀,z₃=f−f₀ とすると

\[
 f_2=2yz_1,
\]
\[
 f_1/2=z_2+z_3+(2d_0+z_1)z_1,
\]
\[
 H=x^2z_2+y^2z_3+(2y^2d_0+2y^3+y^2z_1)z_1.
\]

最後の二式から z₁ 項を除いたものを X,Y とすれば

\[
 z_2=(Y-y^2X)/\Delta,\qquad
 z_3=(x^2X-Y)/\Delta.
\]

従って (f₂,f₁,H) と (z₁,z₂,z₃) は、この既知単元の localization で同じ ideal を生成する。未知の因子を割っていない。

## 最後の一方程式と exact ideal equality

4×4 Newton 恒等式から、全一般行列 C に対して

\[
 f_3+4(\det C-\Delta^2)=\tfrac12f_1(q_1+2q_3)
\]

が成立する。この全係数恒等式も検査した。従って最後の条件は det C=Δ² である。

C₀ を b=e=1,d=d₀,c=c₀,f=f₀ の行列とし、

\[
 \boxed{\mathcal F(x,y,a,g,h,i,j)
      =y^2\bigl(\det C_0-(x^2-y^2)^2\bigr).}
\]

と定める。これは**七変数の77項整数多項式**で、通常の全次数は8である。gauge 固定後の座標では斉次とは限らない。完全な係数は証明書に保持する。graph へ戻すと

\[
 f_1=f_2=0,\qquad f_3=-4\mathcal F/y^2,\quad
 f_4=4\mathcal F/y^2.
\]

加えて det C−Δ² を d,c,f の順に graph へ代入する telescoping 差分を計算し、各差分を z₁,z₂,z₃ の明示的な倍数として保存した。これと前節の逆式により

\[
 (f_1,f_2,f_3,f_4)
 =(d-d_0,c-c_0,f-f_0,\mathcal F)
\]

が QQ[x,y,a,c,d,f,g,h,i,j,1/(xyΔ)] 上で成立する。b=e=1 の gauge を入れた後の等式である。GB、因数分解による不要因子の選別、標本からの同値性推測は使わない。

## chart が実際の像成分と NEW25 に交わること

新しい β 点は取らず、既に使用・独立検証済みの **同じ NEW25 点[0] の C,D**だけを使う。D について

\[
 S=\tfrac12\operatorname{tr}D^2,\quad P=\det D,\quad
 \mathfrak d=S^2-4P
\]

を計算すると P≠0,𝔡≠0。従って D は regular semisimple である。

η²=𝔡 の二次 étale algebra 上で X=(S+η)/2、Y=(S−η)/2 を D² の二固有値とし、

\[
 P_X=(D^2-YI)/(X-Y)
\]

を作る。projector の冪等性と P_X+conjugate(P_X)=I を厳密に確認した。D を対角化したとき、X=x² の block における b,c の積は

\[
 bc=\tfrac12\operatorname{tr}(P_XCP_XC)
    -\frac{\operatorname{tr}(DCP_X)^2}{4X}.
\]

第一項は a²+bc、第二項は a² である。共役側は ef になる。この二つの積を exact quadratic arithmetic で計算すると

\[
 \boxed{(bc)(ef)\in\mathbb Q^\times.}
\]

完全な有理数は証明書に保存した。分子295 bit・分母282 bitで、補助的な32003還元値は10234。従って b,c,e,f はすべて非零であり、有限平方根 gauge が本当に適用できる。この証明は固有ベクトルを近似したり、別の β 点を探したりしない。

同じ点で前 packet の四式の勾配 rank4 も受理済みである。従ってこの chart は実際の16次元像成分の非空開集合に交わり、その正規化した切片の選ばれた成分は6次元である。さらにこの点は NEW25 点なので、NEW25 をこの chart の外に捨ててしまうこともない。

## 有限な軌道切片という範囲

regular D の順序付き対角化の違いは maximal torus と有限 Weyl の選択である。b=e=1 は torus の二自由度を平方根で固定し、残る effective 相対符号は有限である。従って、選ばれた像成分の一般軌道は、この七座標 hypersurface の対応する成分上の**有限個の正規形**で表される。b,c,e,f が非零の一般部では、全 Weyl 選択8と effective sign2で高々16個となる。

これは一般軌道パラメータへの6次元の有限 algebraic cross-section を与えるが、商が有理であること、77項多項式が全ての必要な像成分を尽くすこと、その全域の既約性・正規性は証明していない。

この開集合で、D の中心化 Lie algebra は二次元 torus、b,e≠0 により C をも固定する Lie algebra は0である。従って同じ NEW25 点の orbit は10次元で、一般 NEW25 の quotient は5次元となる。NEW25 の未知条件は、周囲の6次元切片の選ばれた成分上の**さらに一つの未知 divisor**である。今回の 𝓕 自体を NEW25 方程式と呼ばない。

## 保存と再現

* `code/q_firstjet_two_adjoint_normal_form_20261103_v63.py`
* `reports/Q_FIRSTJET_TWO_ADJOINT_NORMAL_FORM_20261103_v63/{report.json,certificate.json.gz}`

全四式、graph、逆の線形式に用いる係数、determinant telescoping cofactors、77項の全係数、同じ点の D² projector と bc,ef の完全な値を保存した。

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_firstjet_two_adjoint_normal_form_20261103_v63.py \
  --output work/two-adjoint-normal-form-replay
```

60秒・2GB上限。16変数の展開、新しい β 点、高次数不変式 probe、GB、NEW25 の新しい方程式は含まない。
