# actual42上のEuler reduction：degree28のscalar targetとcofactor障害

2026-11-03。自由formal45でのconstant-span障害を繰り返さず、実normal42への非線形代入を計算する。
新しいEuler型Red circuitを明示し、そのactual moment determinantの次数が一般beta・generic moment parameter a上で正確に28であることを証明する。
未知残余の空性そのものは未証明である。

## 1. 元の全代数を保つ新しいnormal circuit

既証明の全QQ cross identityを z_k*n_j=P_kj(n) とする。
normal次数dのhomogeneous polynomial f_dにEuler公式を使えば

    z_k*f_d = (1/d)*sum_j P_kj(n)*partial_j(f_d)

が元の全S61 quotientで厳密に成り立つ。
これをhomogeneous成分ごとに繰り返し、z-wordの外側からの順序を昇順と固定する。
Pはlinear+quadraticなのでRed(z^m*f_d)の次数はd以上d+m以下である。
plane cofactor恒等式にこのRedを用いて

    S_l^E=−Red_E(sum_i B_li(z)*K_i(n)),
    H^E=−Red_E(sum_i B_gi(z)*K_i(n))

を定義する。S_l^Eの次数は4以下、H^Eは5以下である。
これは以前の凍結Red実装のbyteやcircuitを変更したものではなく、新版の実42多項式recipeである。
各reduction stepは元cross関係の帰結なので、旧global quintic証明の二つの逆写像はこの版でも成り立つ。
従って全134式、全support、nilpotentを保つ同じactual代数を表す。

## 2. actual次数の正確な上界28

remainder r_iのh^4係数はS1,S2,nに一次で、Hを含まない。
全134についてこれを保存QQ presentationで確認する。
一般betaでも、q_i(h)のh次数が4以下、元式が二次、monic quinticでのdivisionではh^8からH*h^3までしか生じないため同じ構造が成り立つ。

moment r(a)=sum a^i r_iを r=(L0+Q,L1,L2,L3,L4) と書く。
actual normal次数はQ≤8、H≤5、L0,...,L3≤5、L4≤4である。
universal companionから作るbottom4 determinant deltaについて、この重みでの正確な次数上界は28、その最高部分は

    H² L3² L4² − H L3³ Q − L4 Q³

となる。formal weight7の最高項−L4 Q³だけをactual最高項と誤認しない。
この三項全部を厳密integer polynomial determinantで検査する。

## 3. 実nonlinear substitutionによる非零QQ証明

normal方向e0に制限し、S1,S2の四次部分とHの五次部分をfinite Euler circuitで評価する。
K_iの二次部分をf、P_kの二次vector fieldをV_k、JacobianをJ_kと書くと、必要な三次cofactor評価は

    D_a D_b f = Hess(f)(V_a,V_b) + grad(f) J_b V_a,

    D_a D_b D_c f
      = Hess(f)(J_b V_a,V_c) + Hess(f)(V_b,J_c V_a)
        + Hess(f)(V_a,J_c V_b)
        + grad(f) Hess(V_c)(V_a,V_b) + grad(f) J_c J_b V_a

だけで計算できる。四次・五次42変数多項式を展開しない。
Euler分母は二次cofactorで2*3、三次cofactorで2*3*4である。

QQ cross inverse、全K/P係数、plane cofactorのp-integralityと全modular modelへの一致は既に独立監査済みである。
従って上の有限有理数circuit全体はZ_(32003)上に定義される。
今回の実際のreduction値は次の通りである。

    S1_4=15338, S2_4=11727, H5=31338,
    Q8=30427, L4_4=2943, L3_5=13746,
    delta28=31492 != 0             (mod32003).

a=0のr0で非零なので、generic aでもこの最高係数は非零である。
非零像を持つZ_(p)の元はQQでも非零であり、これは**actual42のQQ係数が非零である証明**である。
有限体で零だった式をQQの零式に昇格したり、有限体unitから残余の空性を推論したりしない。
巨大なQQ値を展開保存したとも主張しない。正確なQQ値はpinした有限circuitで一意に指定される。
一般betaへの移行も、非零有理係数の一つのp-unit特殊化によるものである。
これでactual delta_E(n,a)の次数は正確に28となる。

## 4. 実42だけのexact membership target

既知のnormal scalar idealを

    J0=(H3, compatibility2583, cleared-plane210, commutator42)

とする。全2838生成元の次数は4以下で、既証明の等式は

    J0 = m intersect I_B = ker(R42 -> actual A)

である。Euler版でもR42→Aは同じnormal mapなので、このidealは変わらない。
delta_E(0)=0であり、sharp moment証明によりdelta_Eは未知残余B(a)上の非零因子である。
従って正確に

    B=0  iff  delta_E(n,a) belongs J0

となる。これはformal45だけのnegative statementではなく、全actual42上のscalar membershipである。
必要性はB=0ならJ0=m、十分性はdelta_EがBをannihilateし、しかもB上でregularであることによる。

さらにdelta_Eの次数が28、生成元次数が4以下なので、polynomial cofactor表示

    delta_E=sum_j c_j(n,a)*J0_j(n)

が存在するなら、少なくとも一つのc_jのnormal次数は24以上でなければならない。
これは**actual42上の低次数cofactor探索への具体的障害**である。
次数24以上でのmembershipや、より高いdegreeを用いた尽くしを排除していない。
cofactorの上界も与えない。次に必要なのは、この明示scalar targetのQQ certificate、または別の実際のnormal構造である。
