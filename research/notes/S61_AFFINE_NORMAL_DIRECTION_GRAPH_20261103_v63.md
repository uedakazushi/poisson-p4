# 定数commutatorで全compact normal-direction graphを解く

2026-11-03。全134式を保持した未知残余を42個のprojective normal方向だけの式へ移す。compact normal-direction projectionは閉埋込みとなり、各標準chartでは4個のfiber座標を全て多項式で消去できる。アフィン未知残余の空性は未解決である。

## 1. 全sourceとprojective fiber

有限5生成packetの記法を使い、z_k*n=P_k(n)=T_k*n+Q_k(n)、残る三式H=H1+H2、五式F=f_plane(z)+K1(n)+K2(n)とする。P,Q,H,Kの添字1/2はnormal次数である。全134式と定数可逆行操作で同値な表示である。

[w:z:n]と[v]∈P41についてn=r*vと書く。129個のhomogeneous Gは(rを共通因子として除くprojective incidenceでは)

    zeta_k*v=s*T_k*v+Q_k(v),
    s*H1(v)+H2(v)=0,

となる。ここでs=w/r、zeta=z/rである。r=0の点は存在しない。実際w≠0なら共通T-eigenvectorが必要となり、可逆な[T0,T1]に反する。w=0ならcrossrank126からzも零となりprojective点ではない。これは五式を落としてよいという意味ではない。

## 2. 全42chartでの明示的なBezout

C=[T0,T1]は既にQQ p-unit determinantを持つ。v_i=1のchartでlambda_k=(T_k*v)_iと置き、

    U_k=lambda_k*v-T_k*v,
    V_k=Q_ki(v)*v-Q_k(v)

とする。混合compatibilityはs*U_k+V_k=0である。行ベクトル

    u0=row_i(C^-1)*(T1-lambda1 I),
    u1=-row_i(C^-1)*(T0-lambda0 I)

はu0*U0+u1*U1=1を満たす。これは次の一様なhomogeneous恒等式から従う。

    v_i²*C*v=(v_i*T1-lambda1 I)*U0_hom
                  -(v_i*T0-lambda0 I)*U1_hom,
    U_k_hom=lambda_k*v-v_i*T_k*v.

C^-1を掛けi成分を取るとv_i³である。producerは三つのiでこの全三次係数をmod32003で検査し、普遍的な行列展開の証明をここに与える。係数がQQで正しいことは恒等式の形式と保存QQ cross inverseの定義による。modp標本から一般恒等式を推測していない。

従ってB_i=u0*V0+u1*V1とおけば

    s=-B_i(v),
    zeta_k=-B_i(v)*lambda_k+Q_ki(v)

となる。B_iの次数≤4、zetaの次数≤5である。fiber変数を反転する局所化やcyclic determinantは必要ない。三つ目のU/V、三本のH、五本Fは全て残す。

## 3. exact compact graphと実アフィン残余

方向chartのidealは、123本のV_kj−B_i U_kj、三本のH2−B_i H1、五本のhomogeneous Fに上のs,zetaを代入した式で与えられる。次数の上界は順に6、5、10であり、係数は小さな有限逆行列recipeで指定する。巨大QQ多項式の全展開を保存したとは主張しない。

これらの式の商環と元のcompact incidence chartは、明示的な逆写像を持つ。従ってcompact Y→P41は**スキームとして閉埋込み**である。projective incidence全体がbeta基底上でproperだが、実アフィン残余はs≠0、すなわちB_i≠0の開部分であり、これをbeta上properと混同しない。実qの回復は

    n=-v/B_i,   z=lambda-Q_i/B_i

である。全42標準chartを用いるので未知点を失わず、nilpotent構造も保持する。

残った具体的目標はこの方向idealでB_iのnilpotenceを**QQで**証明すること、または真のs飽和を制御することである。modp nilpotenceだけでは一般QQへ移らない。properな族でも点が固定hyperplaneから動く反例があり、flatnessや正確なgraded span等の追加証明が必要である。

## 4. compact Yの空性は目標にできない

凍結toric中心ではqの第8,16,24,32座標の単位軸（zero-based45q）がw=0で元134式の全最高q式をQQで満たす。対応するfree uは11,22,33,45である。元sourceで各軸の唯一の平方項は2u11²u36、2u22²u27、2u33²u13、2u45²u5で、これらbeta係数はこの中心で零である。一般betaで同じ軸が残るとは主張しない。

より強く、既存 `S61_AFFINE_BOUNDARY_DIMENSION_20261102` は全135式の未飽和一般無限遠の次元が**正確に3**と証明している。binary/flag incidenceがbetaを優越するため、betaを変えてもgeneric raw compact Yは非空である。従って別のmodp中心でcompact Y全体の空性を証明してgenericへ移す方針は、このsourceに対して成立しない。無限遠のs=0成分を含むraw Yと、s≠0の未知残余を区別する。

## 5. 保存と再現

`code/s61_affine_normal_direction_graph_20261103_v63.py`、`reports/S61_AFFINE_NORMAL_DIRECTION_GRAPH_20261103_v63/`。先のnative QQ出力のhashを確認して全P係数をmodpへ落とし、全normal係数・H・五式のK・plane Fを小さな圧縮モデルとして保存する。QQの一般recipeは元source、cross inverse、constant commutator inverseで一意に指定される。上の四無限遠点と全三次Bezout係数を検査する。

全S61一般五点尽くし、アフィン残余の空性、QQでのB_i nilpotenceはなお未解決である。
