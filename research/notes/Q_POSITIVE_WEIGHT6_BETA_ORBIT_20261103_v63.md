# 正軸のintrinsic一次ジェット軌道は11次元：T21などの全閉包を排除

2026年11月3日。正weight6軸のprimitive betaを全S61座標から再構成する。
ここでSp4は値のsymplectic frameの変更、Cstarは一次ジェットの一様scalar作用である。
後者は同次座標のdiag(1,t,t,t,t)と全bivectorのt^{-2}倍により得られ、
S61のbetaをt^{-1}倍する。したがって既知成分とその一次像閉包はこの両作用で不変である。

QQ[eta]でSp4 Lie algebraの全10方向を `A*J+J*A^T=0` から別に作り、
intrinsic二形式への有限pullback作用の微分

    −A^T*beta−beta*A−sum_i (A*z)_i * d_i beta

の全24係数を計算する。
beta自身のscalar列を加えた24×11行列の非零11×11 minorは
非零有理定数倍のetaの一つの冪となる。証明書は全行列・全Lie basis・minorの行番号と値を保存する。
従ってeta≠0の全点でSp4×Cstar軌道は正確に11次元である。
群次元も11なので、単なる標本による下界ではなく実際の軌道次元である。

既存の独立監査済み全一次像の次元表から、次の全閉像は11次元より小さい。

|既知成分|一次像閉包の次元|
|---|---:|
|T21|10|
|CA,CB|9|
|zero-log 5H|10|
|C−磁気|10|

もし正軸の一点がこれらのactual成分閉包に入れば、そのbetaの全11次元軌道が
該当する閉一次像に含まれることになり、次元に反する。
従って非零正軸はこれら五成分の全閉包の外にある。
CA/CBと後二つは別packetの独立のobstructionでも排除済みだが、T21排除は本稿の新しい帰結である。

T21はMPSの内部component IDであり、R19の境界記号と混同しない。
T04,...,T38の他の成分の像は次元11以上なので、この議論から排除しない。
新成分・全成分数の増加・全枝の尽くし上げを主張しない。

* `code/q_positive_weight6_beta_orbit_20261103_v63.py`
* `reports/Q_POSITIVE_WEIGHT6_BETA_ORBIT_20261103_v63/{report.json,certificate.json}`
