# 全射影 CA/CB に共通だが N に入らない階数二点：独立監査

2026年11月3日。frozen not-N producer の全 QQ 証明書を別 directory で再実行し、秒数以外一致した。
独立 sparse150 canonical 化と実際の二つの全域 family の基底を使い、さらに **CBへの所属自体** を
Bruhat table に依存せず確認した。

点は P=∂y∧V、V=(0,4x²,0,−2xz,3h²+z²+2xw) である。
CAのg=0表示 f=−4x²,H=3h²z+2wxz+z³/3 がこれを与える。
同じ canonical150 vector は、CBのb=0基底(f3,H16の順)の係数

    (0,6,0, 0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0)

を permutation (1,3,2,0,4) で押し出したものと厳密に一致する。
permutation は old index→new index の規約で、全150係数を比較した。
従ってこれは実際の全射影 CA∩CB の rank≤2 点である。

intrinsic 条件 E∧v∧P=0 の全175×5行列は rank4、核は v=∂y の一線だけ。
canonical bivector の定数 covector contraction map は rank5で核を持たない。
一方 N の f=0境界の全 raw tensor は既に canonical で h,x の二 covector を零化するため、
この点はその境界に入らない。

残る N_g0 では unique factor L を回復して U=V/L 上の canonical quadratic field C を得る。
C は線形 conjugacy と全体倍率を除いて intrinsic であり、N_g0 の表示なら独立な ell,k が

    C(ell)=lambda*ell²,
    C(k)=−alpha*ell²−4lambda*ell*k

を満たす。これは raw Y の divergence=−beta*ell と
C=Y−div(Y)E_U/5 から従う。
独立コードで点の C を再構成し、第一条件で lambda=0 は kernel rank4から不可能、
lambda≠0 は w²,h²,z² の係数により ell∝x、lambda*ell=12x/5と確認した。
第二条件は C(k)+48xk/5 が x²の倍数であることを要求するが、その全係数は k∝x を強制する。
二 covector の独立性と矛盾するため、点は N_g0 にも入らない。

受理済み全N境界=N_f0∪N_g0 により、この点は **全射影 N の外**にある。
従って全域 CA∩CB=N は成り立たない。S61内・rank4開部分での既存の等式とは矛盾しない。
この一つの点だけから例外族の次元や全域交叉の全成分数は決めない。
証明書は code/ambient_rank2_exception_review_20261103_v63.py と同名 REVIEW report に保存する。
