# nil22-image・base括弧零のH+既約quartic全支持は射影次元24以下

2026-11-07。主義務Gの型6、\(N=y\partial_x+w\partial_z\) を扱う。
全27自由係数源でbase括弧 \(f=\{y,w\}=0\) とした**全24係数**を保持する。
Pfaffianの最高四次項を零とする既存の重複H用14次元定理は使わない。
全Jacobi idealの無期限Groebner計算も使わず、平面free divisorの構造から次元を評価する。

**定理（今回証明、独立読解監査対象）。**
この \(f=0\) 型で実退化因子が単純 \(H\) と既約quarticの和となる全点の、
PGL商の前の射影飽和像は次元24以下である。
従って射影次元25,26,27に限られる未知G優越成分の一般元は、この型にはならない。
既知成分へのE全点帰属は主張しない。

有限契約は jobs/hyperplane_quartic_nil22_cotangent_g_bound.json、
実行器は code/classification_quartic_nil22_cotangent.py、
算術証明書は CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_COTANGENT_G_BOUND_20261107.json。
算術と本文の幾何学的証明を区別する。

## 1. 元の全源は二つのquadratic anchorと任意の磁気項

入力は CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json のrow5である。
全110線形行から全90係数・四Jacobi・実Pfaffianを独立Fraction算術で再構成した。
\(f=0\) は正確に \(n_{18}=n_{19}=n_{24}=0\) であり、24自由係数を残す。
fiberを \(q=(x,z)\)、baseを \(p=(y,w)\) とすると

\[
 P=\partial_x\wedge U+\partial_z\wedge V+
       (\alpha x+\beta z+R)\partial_x\wedge\partial_z, \tag{1}
\]

\(U,V\) は \(p\) だけに依存する次数二以下のベクトル場、
\(R\) は任意の次数二以下の式である。
実際

\[
 R=n_3+n_4y+n_5w+n_6y^2+n_7yw+n_8w^2
\]

で、これら六磁気係数は全Jacobiから消える。
divergence条件は

\[
 \alpha=w+\operatorname{div}V,\qquad
 \beta=-y-\operatorname{div}U. \tag{2}
\]

を強制する。\(\Delta=\det(U,V)\) は実affine Pfaffianである。

\[
 K=[U,V]-\alpha U-\beta V
\]

と置くと、全四Jacobiは符号を含めて
\((-K_y,0,K_w,0)\) に等しい。
さらに、全24係数を自由にしたまま

\[
 U(\Delta)+y\Delta=\det(U,K),\qquad
 V(\Delta)+w\Delta=\det(V,K) \tag{3}
\]

が恒等的に成り立つ。これは全Jacobiを課したときのlogarithmic条件を、
標本計算ではなく元の全式のsyzygyとして示す。
全式・磁気自由度・(2),(3)を別の疎多項式実装で展開した（約0.184秒、RSS13,816 KiB）。

## 2. H+quarticは平面free quintic、exponents(2,2)

base座標を \((h,y,w)\) として全式を斉次化する。
\(H=(h=0)\)、実quarticを \(F(h,y,w)\)、平面五次を \(D=hF\) と書く。
仮定により \(F\) は既約四次で \(h\) と異なり、\(D\) は被約である。
\(U,V\) は \(h\) を消す斉次quadratic derivationであり、(3)は

\[
 U(D)=-yD,\qquad V(D)=-wD
\]

を与える。Euler場を \(E=h\partial_h+y\partial_y+w\partial_w\) とすれば

\[
 \det(E,U,V)=h\det(U,V)=hF=D.
\]

従ってSaito criterionにより
\(\operatorname{Der}(-\log D)\) は \(E,U,V\) をbasisとする自由加群で、
projective exponentsは \((2,2)\) である。
これは任意の平面quarticに先に課した追加仮定ではなく、実Poisson条件からの帰結である。

\(\nabla D\) のsyzygyは
\(U+(y/5)E,\ V+(w/5)E\) の二つで自由に生成され、いずれも係数次数2である。
従ってJacobian idealは

\[
 0\to\mathbb C[h,y,w](-6)^2
       \to\mathbb C[h,y,w](-4)^3\to J_D\to0
\]

を持つ。Hilbert numerator \(1-3t^4+2t^6\) の二階微分から、
平面五次の全Tjurina数は正確に

\[
 \tau(D)=12. \tag{4}
\]

## 3. quarticは有理、全Milnor数は5または6

\(F\) の正規化種数を \(g\)、各特異点の枝数を \(r_p\) とし
\(B=\sum_p(r_p-1)\) と置く。
plane quarticの算術種数は3なので

\[
 \mu(F)=2(3-g)-B=6-2g-B.
\]

\(k=\#(H\cap F)\) は**相異なる交点の数**とする。
局所公式
\(\mu(fg)=\mu(f)+\mu(g)+2I(f,g)-1\) とBézoutから

\[
 \mu(D)=\mu(F)+8-k.
\]

\(\tau(D)\le\mu(D)\)、(4)、\(k\ge1\) を合わせると

\[
 2g+B+k\le2. \tag{5}
\]

従って可能性は正確に

| \(g\) | \(B\) | \(k\) | \(\mu(F)\) |
|---|---|---|---|
|0|0|1または2|6|
|0|1|1|5|

となる。この整数リストも全列挙した。

## 4. 必要な小Milnor数の局所型は有限でTjurina=Milnor

ここで全plane singularityの分類定理は不要である。
\(\sum\delta_p=3\)、\(B\le1\) だけから、次の局所normal-form補題で十分である。

* 既約一枝で \(\delta\le3\) なら \(A_2,A_4,A_6,E_6\) のいずれか。
* 二枝で \(\delta\le3\) なら \(A_1,A_3,A_5,D_5\) のいずれか。

一枝の場合、multiplicity2ならWeierstrass表示を平方完成して
\(u^2+v^{2\delta+1}\) となる。
multiplicity3の最小deltaは3であり、等号ではPuiseux pairが(3,4)となる。
主部 \(u^3+v^4\) のJacobian quotientの全monomialのweightは主部より小さいので、
高weight項はJacobian idealに入り、有限決定性と有限回の座標消去により \(E_6\) を得る。
multiplicity4以上なら \(\delta\ge6\) で除外される。

二枝が共に滑らかなら、その接触次数 \(m=1,2,3\) により
\(u^2+v^{2m}\)、すなわち \(A_1,A_3,A_5\) となる。
一枝が特異なら、deltaの加法
\(\delta(C_1\cup C_2)=\delta(C_1)+\delta(C_2)+I(C_1,C_2)\)
から、唯一の残りはordinary cuspと、それに横断する滑らかな枝である。
これは \(u(u^3+v^2)\)、すなわち \(D_5\) となる。
両枝が特異、または三枝以上はdelta/Bの上界に反する。

これらのnormal formは全てweighted homogeneousなので、
Euler関係により函数自身がJacobian idealに入り、\(\tau=\mu\) である。
従って \(F\) の全局所analytic型は有限リストの有限組合せで、

\[
 \tau(F)=\mu(F)\in\{5,6\}. \tag{6}
\]

未確認の高次数cuspidal curveのmoduli分類を引用してはいない。
各型は有限決定的なので、これらの組合せは有限jet条件による有限個のconstructible strataで扱える。

## 5. quarticのequianalytic familyは次元9または8以下

\(\Sigma\) を \(F\) のTjurina zero-dimensional schemeとする。
長さは5または6である。以下で

\[
 H^1(\mathbb P^2,I_\Sigma(4))=0 \tag{7}
\]

を示す。

長さ5以下のzero-schemeは一般のregularity上界 \(\operatorname{reg}I_Z\le\operatorname{length}Z\)
から (7) を満たす。この上界はline residual完全列による長さの帰納法で得られる。
長さ6でschemeが一直線上になければ、
あるline \(\ell\) が \(2\le m=\operatorname{length}(\Sigma\cap\ell)\le5\) を満たす。
異なるsupportを結ぶline、または一点支持ならtangent lineを取ればよい。
残余schemeの長さは \(6-m\le4\) なので
\(H^1(I_{\operatorname{Res}_\ell\Sigma}(3))=0\)。
また \(H^1(\ell,\mathcal O_\ell(4-m))=0\)。
residual完全列が (7) を与える。

最後に長さ6の \(\Sigma\) がline上にあることは、既約quarticでは起こらない。
全三つの偏微分はそのline上で長さ6のschemeを消す三次式なので、
全偏微分がlineを因子に持つ。Euler公式から \(F\) 自身もlineを因子に持ち、矛盾する。

次数4の平面曲線のparameter空間は \(\mathbb P^{14}\) である。
固定した上記有限analytic型のequianalytic stratumの接空間は
\(H^0(I_\Sigma(4))/\mathbb C F\) に入り、(7)により次元は \(14-\tau(F)\) 以下。
特異点を固定したtangent条件ではなく、位置の移動も含むTjurina条件である。
従って全quartic familyは

\[
 \begin{cases}
 \dim\le8,&\mu(F)=6,\\
 \dim\le9,&\mu(F)=5.
 \end{cases} \tag{8}
\]

有限analytic型の和を取っても上界は増えない。

## 6. line Hの選択を戻しても平面pair familyは次元9以下

\(\mu(F)=6\) なら (5) より \(k\le2\)。
固定した既約quarticに対し、一般lineは滑らかな点で四つの相異なる交点を持つ。
従って \(k\le2\) となるlineの集合はdual \(\mathbb P^2\) の真の閉部分集合に入り、
次元は高々1である。(8)からpair \((H,F)\) の次元は高々 \(8+1=9\)。

\(\mu(F)=5\) なら \(k=1\) であり、固定した \(F\) に対するHの選択は有限である。
交点が滑らかなら接触次数4のlineなので、char0のHessianによる有限のinflection lineに入る。
特異点は有限個で、そのmultiplicityは3以下。
従って全交点数4を一つの特異点に集中するlineは、その点の有限なtangent方向だけである。
multiplicity4の点を持つ既約quarticは存在しない
（その点を中心とする次数4のbinary formになり、複素数体上で分解する）。
よってこの枝でもpair family次元は高々9である。

以上で

\[
 \dim\{(H,F)\text{ が今回の全条件を満たす平面pair}\}\le9. \tag{9}
\]

全pairの分類や各stratumのexact dimensionは必要ない。

## 7. 固定した実divisorのPoisson fiberは射影次元9以下

既約quarticを含む \(hF\) の本質変数数は3である。
固定した本質三次元base、固定したpair \((H,F)\) について考える。
Saito basisが \(E,U,V\) で自由なので、任意の斉次quadratic logarithmic場は

\[
 a_1(h,y,w)E+bU+cV
\]

と書ける。\(h\) を消すなら \(a_1h=0\) より \(a_1=0\)。
従ってquadratic anchorの空間は正確に二次元であり、
rank4の二anchorはそのGL2 frame、四parameterで被覆される。

anchorを固定するとJacobiはfiber-fiber括弧の \(x,z\) の二係数を一意に決める。
任意に残るのはbaseの斉次二次磁気式六parameterのみである。
単純因子Hに対するh-Casimir代表を取り、base座標を局所的に一つ選んで数えており、
同じdivisorに対する座標の選択を新たな幾何学的parameterとして加えない。
fiber座標のbase一次式による加算は追加parameterを要しない。
実際 \(x'=x+\ell,\ z'=z+m\) とすると磁気項は

\[
 R'=R-\alpha\ell-\beta m+U(m)-V(\ell)
\]

へ移り、これも既に全六次元の磁気二次式の中にある。
GL2 fiber変更はanchor frameの四parameterに含まれる。
従って固定projective divisor fiberは高々 \(4+6\) 次元の係数族で被覆され、
非零tensor scaleを除いて高々9次元となる。

## 8. 全P4飽和像とGへの帰結

本質三変数のlinear subspaceの選択は
\(\operatorname{Gr}(3,5)\) の六次元である。
その中の平面pair familyは (9) により高々9次元。
従って全divisor familyは高々15次元である。
第7節の固定divisor fiber上界を足して

\[
 \boxed{\dim(\text{型6 }f=0,\ H+\text{既約quartic全射影支持})\le
          6+9+9=24}. \tag{10}
\]

これは全支持の上界であり、generic source、最高部だけの条件、
一つのnormal form又は有限標本からの推論ではない。
未知G優越成分の次元25,26,27という既済の命題から、この型はその一般元にはならない。

非零base括弧の全支持は TWO_QUADRICS_NIL22_BASE_SUPPORT.md により、
重複Hまたは既知23次元H_QUARTIC_ESSENTIAL5へ処理済みである。
従って今回の結果を合わせると、型6の全H+既約quartic支持がGの一般型から除かれる。
f=0の全点の既知成分帰属、新成分の有無、Eの完成は別の未解決問題として残す。
