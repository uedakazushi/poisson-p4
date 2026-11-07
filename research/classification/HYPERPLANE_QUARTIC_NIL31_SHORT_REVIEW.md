# H+quartic・nil31短translation全支持の次元上界：独立検証

2026-11-07。対象は型4の全18parameter source、
\(Z=y\partial_x+z\partial_y+\partial_w\) である。
根元の全90係数からJacobi・Pfaffianを再構成し、
code/classification_quartic_nil31_short_review.py でproducer/CASをimportせず検証した。
現行証明書は reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_REVIEW_V2_20261107.json、
契約は jobs/hyperplane_quartic_nil31_short_review.json。
実行約0.216秒、最大RSS13,952 KiB。gateは入力契約だけの検査であり、
以下の分岐・次元論は独立した数学的確認である。
初回のimmutable証明書は保存し、V2では同じ九枝・全行列を持つ公開producer report
CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_SUPPORT_20261107.jsonへ入力を移した。
元のprivate work pathは現行jobの依存ではない。検証算術とGroebner不要の証明は同じである。

## 1. 全sourceとcoreの対応

既存の独立Fraction算術を用い、110個のdivergence/Lie線形行、
全90係数、六括弧、四つの全Jacobiと実Pfaffianを別に再構成した。
元sourceの自由度18、Jacobi係数をmonic重複除去した全13式を確認した。
coreは

\[
 j=n_{10},\quad a=n_{11},\quad k=n_9,\quad
 t=n_{15},\quad b=n_{16},\quad c=n_{17}
\]

だけに依存する次の六式である。これらと元の六core式のmonic集合が一致することも検査した。

\[
\begin{split}
 &(a+b+\tfrac16)t+ck=0,\\
 &-jt-t^2+(b+\tfrac13)k=0,\\
 &c(a-b+\tfrac13)=0,\\
 &cj-\tfrac12ab+\tfrac16a-\tfrac12ct-\tfrac12b^2+\tfrac16b=0,\\
 &\tfrac19j-\tfrac23(a+b)t+\tfrac19t+ck=0,\\
 &\tfrac1{12}j-\tfrac14(a+b)t+\tfrac18t+ck=0.
\end{split} \tag{1}
\]

以下の解法はGroebner基底の完全性を仮定せず、直接全分岐を尽くす。

## 2. c非零の五つの完全枝

\(\tau=ct,\ J=cj,\ K=c^2k\) と置く。
最初、第三、第四式から

\[
 a=b-\tfrac13,\quad
 J=\tfrac12\tau+F(b),\quad
 K=-(2b-\tfrac16)\tau,\qquad
 F(b)=(b-\tfrac13)(b-\tfrac16).
\]

を得る。全六式に代入した残りは、正確に

\[
 F(b)-30(b-\tfrac16)\tau=0,\qquad
 \tau(\tau+2b^2)=0. \tag{2}
\]

一つ目は第五・第六式にそれぞれ \(9c,12c\) を掛けた式、
二つ目は第二式に \(-2c^2/3\) を掛けた式であり、
全六式の代入恒等式をLaurent多項式で確認した。

\(\tau=0\) なら \(b=1/3,1/6\)。
\(\tau=-2b^2\) なら

\[
 F(b)+60(b-\tfrac16)b^2
   =60(b-\tfrac16)(b+\tfrac1{12})(b-\tfrac1{15})=0.
\]

従って \((b,\tau)\) は

\[
 (\tfrac13,0),\quad(\tfrac16,0),\quad
 (\tfrac16,-\tfrac1{18}),\quad
 (-\tfrac1{12},-\tfrac1{72}),\quad
 (\tfrac1{15},-\tfrac2{225})
\]

の五つである。各枝の唯一のcore parameterは \(c\ne0\)。
\(b-1/6\) などで割っていないので、特殊固有値を失っていない。

## 3. c零の全四枝

\(c=0\) では、(1)の第五式の9倍と第一式の6倍を足すと
\(j=-2t\) を得る。残りは

\[
 (a+b+\tfrac16)t=0,\quad
 t^2+(b+\tfrac13)k=0,\quad
 (a+b)(b-\tfrac13)=0. \tag{3}
\]

\(t\ne0\) なら
\(a=-1/2,b=1/3,k=-3t^2/2\)。
\(t=0,k\ne0\) なら
\(a=1/3,b=-1/3\)。
\(t=k=0\) なら \(b=1/3\) または \(a=-b\)。
従って四つの一parameter枝となる。
最後の二枝の交わりを含み、\(c=t=k=0\) の補集合も保持した。

## 4. 全Jacobiの二段線形条件と全rank降下点

残るparameterをlower六変数
\((n_3,n_4,n_8,n_{12},n_{13},n_{14})\) と
magnetic六変数 \((n_0,n_1,n_2,n_5,n_6,n_7)\) に分ける。
元の全13式を各完全core枝に代入すると、core以外は
lowerについての三つの線形条件と、magneticについての四つの線形条件になる。
後者の係数はlower変数に依存しない。
定数項やlowerの二次項は捨てず、dimensionの上界に不要な整合条件として残す。

各行列は元Jacobiから別に抽出し、producer行列と全成分一致を確認した。
Fractionによる有理函数体上のrankだけでなく、
**全最大minorの分子のgcd**を計算し、その全根を分解した。
報告された有理根を除いたgcdに正次数因子が残らないことを検査するので、
未列挙の代数的例外を隠していない。
全例外で両行列のrankを有理数として再計算した。

全九枝のgeneric rankはlowerが3、magneticが4である。
従ってparameterを含むgeneric source次元は

\[
 1+(6-3)+(6-4)=6.
\]

上の五つの \(c\ne0\) 枝、\(c=0,t\ne0\) 枝、
\(c=t=0,k\ne0\) 枝では、許される領域内にrank降下点はない。
最後の二枝の全例外は次の通りである。

| core枝 | parameter値 | lower rank | magnetic rank | 固定parameterの次元上界 |
|---|---|---|---|---|
| \(b=1/3\) | \(a=-1/3\) |2|3|7|
| 同上 | \(a=-1/4\) |3|3|6|
| \(a=-b\) | \(b=1/3\) |2|3|7|
| 同上 | \(b=2/3\) |2|4|6|
| 同上 | \(b=1,1/2,1/4\) |3|3|6|

例外値は零次元なので、この行にはcore parameterの次元1を足さない。
lower方程式が不整合ならfiberは空であり、整合時も次元は \(6-r_L\) 以下。
その上のmagnetic fiberは \(6-r_M\) 以下なので、整合条件の有無にかかわらず

\[
 \boxed{\dim(\text{固定 }Z\text{ の全source})\le7}. \tag{4}
\]

## 5. cofactor軌道と射影全stratum

全20個のinfinitesimal affine変換から、
\(\delta A=[B,A],\ \delta c=Bc-Av\) の全20座標を別に作った。
全20×20行列のFraction rankは13であり、affine cofactor軌道次元は13である。
特定のminorだけを全軌道rankの証明として用いていない。

固定hyperplaneのaffine chartで、(4)とcofactor軌道を合わせた次元は高々20。
hyperplaneの選択は4次元である。
さらに非零scalar \(sZ\) は
\(\operatorname{diag}(s^2,s,1,s)\) によって \(Z\) とaffine共役である。
従ってこの全族はbivectorのscalar倍についてconicalであり、非零部分の射影化で1次元減る。

\[
 \boxed{\dim(\text{型4の全 }H+\text{既約quartic stratum})\le
          7+13+4-1=23}. \tag{5}
\]

複数のhyperplane表示があっても像の次元は増えず、この上界には問題ない。
未知G優越成分の射影次元25,26,27という既存定理と合わせ、
その一般元がこの型になることはない。

これは型4の全点のMPS/T40/他成分への帰属ではない。
また23次元stratumから16次元primitive底への写像が非優越とも結論しない。
新成分数の追加は行わない。
