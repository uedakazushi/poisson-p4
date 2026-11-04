# Finite Tor63像のall-order typed operator：16-generator変換と27 actual C0 maps

2026-11-04。これはsource jetや係数展開でなく、有限unit系と証明済みregular quotientによる正確なoperator構成である。
有限像の零・非零、socleやDのdepthをこの構成だけで決定しない。

## 1. 九つのactual105 unit系

元のpure CA48/CB48とCAのF9を共通physical R0上で使う。
105未知数/105式のprototypeは既存joint105とnode順序まで同じである。
zero7のretained prefix長0,...,6を持つ七系と、二軸(r,0)、(0,s)の系を作る。
prefix長7は元のsource functionそのものを使う。target以外の全演算は元source DAGのhomomorphic cloneである。
各Jacobianは同じ原点105 unit blockなので、新しいGaussianやNewton solveは必要ない。

## 2. 九十個のregular quotientを意味付きで保持する

CA f63とCBの全八zero odd functionの九函数に対し、七coordinate prefix差を対応zero座標で割る。
同座標を零にすると前後のunit系が一致するので、差はその座標に正確に可除である。
63商の原点値はcached first coefficientを用いる。

各plane函数はorder2以上。二軸の値を用いた

```
B(r,0)/r²,
(B(r,s)-B(r,0)-B(0,s))/(rs),
B(0,s)/s²
```

の27商を加える。軸上のorder2と、混合numeratorが両軸で消えることからregularityが従う。
これらは新しいtyped operator qであり、numerator=denominator*qという正確なidentityとregularity根拠を保存する。
商をその分母の零上へfreshに代入し直すことは一度も行わない。
全originはcached source coefficientsで固定し、旧r/gを使ったり評価したりしない。

## 3. 全16-generator unit変換

h=(z7,r²,rs,s²)、G=(f8,g8)と書く。telescopesと二次Hadamardは正確にG=A*hを与える。
CA first7とCB rows18,60,67をselected10に取り、後三本のquadratic matrix Uは原点det625/768。
adjugate3でそのunit inverseを作り、h=T*Gを指定する。
六つのredundant indexに対してq_alpha=e_alpha-A_alpha*Tと置けば、全q_alphaはactual first syzygyである。
eta_u、eta_vもTのtransposeで(-s,r,0)、(0,-s,r)を持ち上げる。

十canonical行Tと六q行を持つ16x16変換のinverseは、最初十columnがA、後六columnがredundant unit vectorsとなる。
一般のU9係数・selected H21係数・redundant A60係数を独立symbolにした整数多項式計算で、
denominator det(U)を払った左右inverseの全512 entryを直接照合する。
これは巨大source展開でなく、有限block inverse identityの検証である。
六cycleの原点も既存preflightの値へ全係数照合する。

## 4. Actual C0 free9への二十七map

CA part of prefix105は、従来のpureCA57 prefixと同じtargetを持つ。
pointed unit解一意性により、plane f63は既存gA、七coordinate商は既存f63 Hadamard係数と一致する。
従ってC0 basis(1,z7*tauA,gA*tauA)へのf_i*tauAの全vectorがそのまま定まる。

cycle q,pの二十七指定pairに対して

```
h_i,beta=q_Ai*p_Bbeta-p_Ai*q_Bbeta,
even=nuA*sum h_i,beta*f_i*f_beta,
odd C0 vector=sum h_i,beta*(aA_beta-aB_beta)*image(f_i*tauA)
```

という九entry operator columnを保存する。八physical rowは両族で同じ順序であり、
共通target内のg_beta*tauB=(aA_beta-aB_beta)+f_beta*tauAを使用する。
これらはraw Tor2の自然像なので全てm_R0²-annihilatedである。

全regularity、unit clone、telescoping、3x3 inverse、Koszul cycle identityは有限operator identitiesで指定・検証される。
generic source functionの多変数係数を展開して零判定したとは主張しない。
imageがD0の640関係に入るかどうか、有限像が非零かどうかは別のmembership gateである。

## 5. 保存と停止範囲

prefix構成checkpoint、全typed scalar DAG、九unit105系、90regular quotient根拠、左右変換matrix、
八cycle、全27 actual C0 image columnsを保存する。serialized size32MiBを実際にguardする。
source5/higher jets、full640 second matrix、oldr/g evaluation、new Gaussian、GBやnew pointは行わない。
二次negative routeは既にinconclusiveで停止済みなので、このoperator構成から高次jetへ自動移行しない。

Code `code/ambient_B18_finite_Tor63_operator_20261104_v63.py`。
Report/folder `AMBIENT_B18_FINITE_TOR63_OPERATOR_20261104_v63`。
