# rank零plane・normal半単純 mixed 三parameter全源の超平面因子

2026-11-08、PRIMARY E。純幾何 producer proof。新しいoperator／kernel／minor／GB／profileは無い。

**定理。** actual被約既約quinticを持つ rank 四 Poisson tensor の rank零plane m二sourceについて、
normal半単純 \(\mathrm{SS\_MIX}\) 三parameter全源は空である。
\(t,s,k\) の全値、全 \(P_1\) 三十方向、全 inhomogeneous \(P_2\) 九方向を保持する。
旧二十三plane／二百三十二curve等の hierarchyを新しい点計算の予算として使わず、
元 first-weight Jacobiだけから actual hyperplane因子を強制する。

## 1. 全 source と二つの bad normal rows

`PLANE_RANK_ZERO_GLOBAL_VECTOR_SEMISIMPLE_FIRST_JETS.md` とその全sourceを受理した
`PLANE_RANK_ZERO_M2_COMPLETION_SEPARATE_REVIEW.md` を入力とする。
全 normalizer と逆を戻した対象の一つの全 sourceは
\[
 P_0=T\wedge E_n+V\wedge(z\partial_z-w\partial_w)+kzw\partial_z\wedge\partial_w,
\]
\[
 H=\operatorname{diag}(1,1,-2),\quad N=b\partial_a,
 \quad T=tH+sN,\quad V=H+N,\quad t\ne s,
 \qquad E_n=z\partial_z+w\partial_w.                            \tag{1}
\]
全 tensorは \(P=P_0+P_1+P_2\)。既受理の canonical全三十源を
\[
 P_1=zL+wM+z^2\mathbf D\wedge\partial_w
  +zw(\mathbf B\wedge\partial_z+\mathbf E\wedge\partial_w)
  +w^2\mathbf C\wedge\partial_z                                \tag{2}
\]
と書く。L,Mはbase-linear base-base bivectors、太字は全constant base vectors。
\(P_2\) は任意の binary quadratic をbase三pairsへ掛けた全九源。
特に \(P_2\) は z-row／w-rowを持たない。

\(\mathbf C=0\) なら full \(P\) の z-rowは全てzで割れ、\((z=0)\) はPoisson不変。
\(\mathbf D=0\) なら同様に \((w=0)\) が不変。
その三次元 hyperplane上でprojective rankは二以下なので、actual quintic Pfaffian はそれぞれz又はwで割れる。
Pf零ならrank四でなく、非零ならactual被約既約quinticに反する。
従って本対象の必要条件は
\[
 \mathbf C\ne0\quad\hbox{かつ}\quad\mathbf D\ne0.                \tag{3}
\]
これは高次 \(P_2\) や他の \(P_1\) を零と置いた議論ではない。
Euler gaugeの weight一六方向も bad rowsを変えないので、(2) の canonical quotientに依存した切断ではない。

## 2. 全 normalSS symbol に共通する first-weight 恒等式

w一 chartで \(u=z/w\)、\(x=(a,b,c)/w\) とし、parameter tと区別する。
\(U=T-V\)、\(R=2V+kI\) と置けば
\[
 \{x_i,x_j\}_0=x_i(Ux)_j-x_j(Ux)_i,
 \qquad\{x_i,u\}_0=u(Rx)_i.                                  \tag{4}
\]
\(u=0\) で \(P_1\) の mixed anchorは \(\mathbf C\)、vertical部分はxにlinear。
元 Jacobi \(\{x_i,\{x_j,u\}\}+\mathrm{cyclic}\) の first-weight、xにlinearな係数を取ると
\[
 C_i(Rx)_j-C_j(Rx)_i-\sum_l C_l\partial_l
       [x_i(Ux)_j-x_j(Ux)_i]=0.
\]
従って全xについて
\[
 \boxed{\ \mathbf C\wedge(-T+3V+kI)x-x\wedge(T-V)\mathbf C=0\ }. \tag{5}
\]
導出では \(P_1\) の他の mixed coefficients の derivative は二つの\(P_1\)の積なので次のweightに属し、
\(P_2\) はさらに低いfiber degreeで、このlinear-in-x係数へ寄与しない。
(4)は projective affine bracketそのものであり、homogeneous raw Jacobiをstrictと仮定していない。
係数体の変数をPoisson中心とした計算でもない。

normal swap \(z\leftrightarrow w\) は
\((T,V,k,\mathbf C,\mathbf D)\mapsto(T,-V,-k,\mathbf D,\mathbf C)\)。
したがって \(\mathbf D\) には (5) の \(V\mapsto-V,k\mapsto-k\) が必要。
二chartと逆を同時に扱い、normal∞での条件を捨てない。

## 3. 一つの nonzero bad vector の全 linear-algebra 条件

\(v=\mathbf C\ne0\)、\(M=-T+3V+kI\)、\(U=T-V\) と書く。
(5) を \(\mathbb C^3/\mathbb Cv\) の exterior squareへ写すと
\(\bar x\wedge\overline{Uv}=0\) for all \(\bar x\)。商はdimension二なので
\(Uv=\beta v\)。元の式は
\[
 \operatorname{im}(M+\beta I)\subset\mathbb Cv.                 \tag{6}
\]
(6) をvに適用すれば \(Vv=\nu v\)、従って \(Tv=(\nu+\beta)v\)。
rank一以下の \(M+\beta I\) のtraceはv上の固有値に等しいので
\[
 3(k+\beta)=2\nu+k,
 \qquad \beta=(2\nu-2k)/3.                                   \tag{7}
\]
さらに商space上で
\[
 \bar T=3\bar V+(k+\beta)I.                                  \tag{8}
\]
全て必要条件であり、未知source係数を反転していない。
\(v=0\) 補集合は既に第1節でfactorとして処理した。

## 4. 全三parameterの badC 支持は二つの line

\(V=H+N\) の固有spaceは、固有値一の \(\mathbb C\partial_a\) と固有値−二の \(\mathbb C\partial_c\) の二つだけ。
二次元Jordan blockのb方向を固有vectorと誤認しない。
従って (7)–(8) は次の全二場合になる。

* \(v\in\mathbb C\partial_a\)：\(\nu=1\)、\(\beta=t-1\)。商spaceのb,c上で
  \(\bar V=\operatorname{diag}(1,-2)\)、\(\bar T=\operatorname{diag}(t,-2t)\)。
  (8) の二固有値を比較すると \(t=3,k=-2\)。sには追加開条件を課さない。
* \(v\in\mathbb C\partial_c\)：\(\nu=-2\)、\(\beta=2-2t\)。商spaceのa,b上で
  \(\bar V=I+N\)、\(\bar T=tI+sN\)。
  (8) のnilpotent係数から \(s=3\)、scalar係数／(7)から \(k=3t-5\)。

つまり全 nonzero \(\mathbf C\) の必要支持は正確に次の二lineの和に含まれる。
\[
 L_{Ca}=(t=3,k=-2),\qquad L_{Cc}=(s=3,k=3t-5).                 \tag{9}
\]
これらの全係数値を含み、generic line又は有限標本だけを見た結論ではない。
\(t\ne s\) のm二domainを戻すと各lineの対応する一点を除くが、以下の矛盾にはその除外さえ不要。

## 5. badD 支持との全交差が空

第2節のactual normal swapを (7)–(8) へ適用すると、全 nonzero \(\mathbf D\) の支持は
\[
 L_{Da}=(t=-3,k=-2),\qquad L_{Dc}=(s=-3,k=-3t-5)               \tag{10}
\]
に含まれる。
全四通りの交差はliteralに空である。

|C支持／D支持|矛盾|
|---|---|
|\(L_{Ca}\cap L_{Da}\)|\(t=3=-3\)|
|\(L_{Ca}\cap L_{Dc}\)|\(t=3\) で \(k=-2\) と \(k=-14\)|
|\(L_{Cc}\cap L_{Da}\)|\(t=-3\) で \(k=-14\) と \(k=-2\)|
|\(L_{Cc}\cap L_{Dc}\)|\(s=3=-3\)|

したがって全 \((t,s,k)\) で \(\mathbf C=0\) 又は \(\mathbf D=0\)。
必要条件 (3) に反し、actual reduced irreduciblequinticの全rank四sourceは空となる。
他のweight一coefficients又は全inhomogeneous \(P_2\) がこのfactorを解消する可能性は無い。
該当normal rowを持たないことを第1節で全sourceのまま確認したからである。

## 6. 全sourceへ戻す矢印と範囲

旧 semisimple first-jet theorem の実normalizer／tensor尺度／zero補集合を戻す。
mixed centralizerでVのsemisimple成分とnilpotent成分が共に非零の全sourceは(1)へ移り、
元の全 \(P_1,P_2\) はinvertible graded変換で保持される。
actual hyperplane因子は実PGL逆でもhyperplane因子なので、対象は全て空である。

VがpureH又はpureN、regularnilcentralizerの全二parameterは別凍結note
`PLANE_NORMAL_SS_TWO_PARAMETER_WHOLE_EXCLUSION.md` で、全saved kernelsとactualshearにより処理した。
本稿の新しいfirst-weight恒等式はそれらにも適用できるが、旧sourceを再計算する必要は無い。
Cartanは既知MPS、common-kernel／low21、normalnil B零等は担当分離した別sourceで、ここから未証明の全pointwise帰属を推論しない。
全rank零plane又は全 E／総component数は未主張。
新component増分は零。旧code／reports／主稿／global／Gitを変更せず、同名manifestへexact inputsを固定する。
