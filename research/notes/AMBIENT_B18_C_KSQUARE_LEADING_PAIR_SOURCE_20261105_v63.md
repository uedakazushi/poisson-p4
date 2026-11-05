# Kbar 二乗の最低次数だけを用いる十分条件

ここでは積の signed 値を計算しない。既存の canonical section による U,S の free44 cycle は原点零である。第一 jet はそれぞれ `L_U v42+d_U r e43`, `L_S v42+d_S r e43` と書ける。v42 は direct section の既存原点列であり、その最後の成分は零、m0 v42=m0 e43=0 である。未知 Gamma2、section1 は Gamma0=0 のためこの第一 jet に入らない。

Product44 の even coefficient は proper J に属する Qeven と a36 の積なので二次以上、odd scalar B は beta7 の一次と b7 の一次の積なので二次以上である。従って二つの cycle の積は free44 で三次以上となる。その三次部分では nuA、gA、入力二次 jet は現れない。`q=Qeven_1 v42_a` と `p=sum z_i v42_bi` とおくと、第一三十六成分は `L L' (q v42_a+p T_A v42_a)`、最後の八成分は `2 L L' q v42_b` と e43 方向の既知項である。T_A は trace mask 0→1,2→3 のみの移動で、nu0=0 を用いる。最後の e43 項は m0e43=0 により今回の次数二抽出から消える。

正確な H を全多項式 m0*c3 に作用させる。sector の n は (1,1,2)、Hadamard は original z0,...,z6 の順序である。v7 の最初の z を一個除き Hi とし、wi=vi+r^n Hi、pair_ij=(Euler_z+2)^(-1)(partial_i wj-partial_j wi) を使う。係数を H の外へ出す S 線形性は使わない。以下の tensor では、R0 単項式と p,q を掛けた後にこの式を適用する。q の八係数は独立な QQ 係数記号として保持する。

積の raw128 representative の次数二は `(0,H(m0*c3)_2)` である。m1*c3 と入力 c4 は H 後でも三次以上、Pi1 と raw 次数二の積も三次以上なので、次数二 class は既存 Pi0 と正確な gr2 projection3 だけで得られる。未知 Gamma2/3、m1、Pi1、lower relation の高次 carry はこの leading class へ漏れない。これは full product class を計算したという意味ではない。

三つの leading class が、old F2 image の二次元平面を張るなら十分である。実際、同平面は `i(v),i(q2)` で張られ、r による gr2→gr3 は rank 一で、`r*i(q2)=i(g)` が gr3 に非零である。従って二つの独立な leading product と、それらの r 倍のうち一つは実 Tnil でも独立である。これは filtration に対する三角性によるので、old→new の strictness は仮定しない。Kbar 二乗像は H0/ellH0 に入り、その次元は三であるから、この十分条件が実値で成立すれば H0=D*K² となる。まだ rank、H0、nil 指数を結論しない。

SOURCE は L_U,L_S を代入せず、q の八係数も対にせず、36 個の二次単項式ごとの線形抽出 tensor を固定するだけである。次の実値 worker に q が必要なら、既存 actual proper first cache と元 radical polynomial source に対する明示的な callback と予算を先に閉じる。
