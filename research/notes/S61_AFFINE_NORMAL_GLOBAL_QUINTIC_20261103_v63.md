# 全normal base上のglobal monogenic quinticと5×670表示

2026-11-03。全S61代数Aは、42-normal base R上で**一つのglobal monic quinticの商**として表せる。以前のquintic-openや1031-generator cyclic-openとはbaseと範囲が異なり、未知qに関するdeterminantの局所化はない。全134式を保持し、実際のR-moduleは5×670の有限表示を持つ。自由rank5やbeta-only有限性、残余の空性は主張しない。

## 1. 小さなplane cofactor証明書

実既知五点のplane代数E5でh=z0はprimitiveである。保存lex表示のz1,z2のquartic補間をp1(h),p2(h)、monic quinticをg(h)とする。元のplane二次式f_i(z)五本に対して

    z1-p1(h)=sum B_1i(z)f_i(z),
    z2-p2(h)=sum B_2i(z)f_i(z),
    g(h)=sum B_gi(z)f_i(z)

の全QQ cofactorを保存した。最初の二つはdegree4 Macaulay行列50×35/rank30、最後はdegree5行列100×56/rank51から得た。cofactor次数は2,2,3、非零項数は31,31,52であり、全三つの多項式積恒等式と分母32003-unitを再検査する。

一般betaでは普遍K-rational plane frameを用いる。generic T01 rootを選ぶ必要はない。frameのoriginは保存三つのfree座標を持つK-rational pointとして選び、凍結中心へ合わせられる。hがE5を分離する条件と保存小minorはbeta-only非空開である。全leading plane quadratic/cubic spanによりcofactorのdegree4/5 boundも保たれるので、同じ有限逆行列recipeを一般betaで使える。

## 2. normal polynomial corrections

全129 crossの既証明identity z_k*n_j=P_kj(n)を使う。全nを含む単項式のz個数を一つずつ減らす決定的なnormal reduction Redを取る。元五本F_iはf_i(z)+K_i(n)へ変わり、K_iは定数項を持たず次数≤2である。

上のcofactor恒等式を全代数へ持ち上げると

    z1=p1(h)+S1(n), z2=p2(h)+S2(n), g(h)=H(n),
    S_l=-Red(sum_i B_li*K_i), H=-Red(sum_i B_gi*K_i)

を得る。次数の上界はS1,S2≤4、H≤5である。これらはR内の実際の多項式であり、新しい独立parameterではない。巨大なnormal展開は作らず、保存QQ cofactorと126×126 cross inverseによる有限arithmetic recipeで指定する。

## 3. 全134のexact univariate reductions

h,z1,z2,nを上式で元の全45qへ戻し、元134式を代入する。各式のh次数は≤8で、monic G(h)=g(h)-H(n)により割ると次数≤4のr_iを得る。producerはS1,S2,Hを**formal placeholders**として扱う自由係数環でも、全134について

    F_i(substituted q)=quotient_i*G+r_i

をQQ係数で全項検査し、全remainderとquotientを保存する。最後にplaceholderを上の実normal polynomial recipeへbase changeする。formal placeholdersを任意の別解族と誤認しない。

R[h]/(G,r_0,...,r_133)とAの間で、h↔z0、z1/z2↔p_l(h)+S_l(n)、n↔nを使う二つのmapは互いに逆である。片側は元134式の代入、もう片側は上の正確なcofactor identitiesによる。したがってこれは全supportとnilpotentを保つ**global algebra isomorphism**である。

C=R[h]/Gはfree rank5だが、A=C/(r_i)はその任意の商である。R-moduleとしてideal(r_i)はh^j r_i (j=0..4)で生成されるから、companion matrix T_hにより5×670の有限表示を得る。このmodule closureはmonic一変数環の恒等式で保証され、複数の非可換な作用行列から自由性を推測していない。

## 4. 残余での純normal多項式回復

別のnormal commutator packetのD84は実際のknown-five idempotent eである。1−D84は定数項が零なので、first-normal factoringでsum c_i(n)n_iと書ける。従って

    z_k(1-e)=sum_i c_i(n)P_ki(n)

で、右辺の次数は≤85である。未知残余Bでは1−e=1なので、全zは明示された純normal多項式で回復する。巨大determinantの展開は不要である。これはnormal閉埋込みの具体的な有限recipeであり、Bが空であることを意味しない。

## 保存と残件

`code/s61_affine_normal_global_quintic_20261103_v63.py`、cofactor `code/data/s61_affine_normal_global_quintic_20261103_v63.json`、receipt `reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/`。全134 h-remainderは小さいplaceholder coefficient modelとして保存する。generic expanded beta/normal coefficientsをserializeしたとの主張ではない。

この新しい5×670 global module、direction-only graph、normal低次数idealは、未知supportを調べる別の有限入力を与える。unit、general length5、全w飽和はいずれも未完了である。
