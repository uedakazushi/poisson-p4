# 有理orbit表示の群座標を平方根なしで戻す

2026-11-04。対象は既受理の第一jet像B_NEW15のみである。五つのsource座標を戻す
有限recipeと、generic Sp4 stabilizerがidentityであることから、残る群座標も
行列式・四則演算で回収できる。全25次元Poisson成分や境界を分類するものではない。

## 有限pullbackと転置の規約

元の24係数は x_k dx_i∧dx_j の係数である。有限作用をgによるpullbackとする。
これはanti-actionで、Lie微分のweightも三つのindexのweightの和となる。
二つの共変量の変換は

    C(g*beta)=g^t C(beta) g^(-t), D(g*beta)=g^t D(beta) g^(-t)

である。spinor行列側をgそのものにしてはならない。元の全Lie係数に対する
転置規約は既受理。今回、非対角symplectic transvectionの有限pullbackを
全24×16 basisと全4352 quadratic coefficientsから作り直し、この向きを検証する。
これは元の点や高次数不変式を評価する計算ではない。

## 群座標の有限有理recipe

入力betaから既受理のrecipeでz=(xi,r,s,q,B)を回収し、b=beta_bar(z)と置く。
同じgeneric orbit上なので、幾何的には一意なg∈Sp4によりbeta=g*bと書ける。
この存在は有理なgを仮定せず、受理済みのproduct mapの支配性とgeneric degree1
またはRosenlicht generic orbitを使う。

16個の未知数T_ijに対して32本の線形式

    T C(b)=C(beta) T, T D(b)=D(beta) T

を作る。共通commutantがscalarなので、解空間は正確に一次元で、T=t g^tである。
非零な15×15 minorを持つchartを一つ選び、残る一座標を1とするCramer式でTを
得る。有限個のchartを辞書式に調べることも完全に指定された有限recipeであり、
未知の多項式恒等式や追加点の探索を必要としない。得たTは可逆である。

標準Jの01成分を用いて

    lambda=(T^t J T)_01

とする。scalar倍symplectic行列なのでT^t J T=lambda J、lambda=t²≠0。
次に元24係数へのfinite pullbackを使い

    v=(T^t)*b

を計算する。三つのcovariant indexを持つためv=t³ betaである。betaの非零な
一座標iを選びk=v_i/beta_iとするとk=t³≠0。したがって

    g=(lambda/k) T^t

である。係数は全て入力betaの有理関数で、平方根は要らない。
g^t J g=J、g*b=betaもそれぞれtの冪から文字どおり従う。
Tの正規化を任意の非零scalar aで変えても、lambda→a²lambda、k→a³kなので
最終gは不変である。

## 十個の有理な群座標

g+Iが可逆なdense openでは

    A=(g-I)(g+I)^(-1), g=(I-A)^(-1)(I+A)

とする。g^t J g=JからA^t J+JA=0が従い、逆も同様である。
A=J S（S symmetric 4×4）と書けるのでSの十個の上三角係数は独立な有理群座標。
元の五つのzと合わせた十五個の有限有理逆座標を得た。
このCayley chartはg=Iで定義され、全Sp4がこの一つのchartに入るとは主張しない。

報告は有限非対角pullbackの向き、scalar指数、Cayleyの全行列恒等式を検証する。
巨大なgeneric intertwinerや15変数逆式の係数展開は行わない。
