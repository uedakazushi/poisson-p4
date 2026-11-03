# 最終次数補足：vector7、square14、135scalar30

2026-11-03。前のbounds8/16/34は有効なまま保持し、companion staircaseで改善する。これはformal45 parameter次数のboundであり、非線形S(n),H(n)への代入後のdegreeを述べない。actual delta、unit、membershipは未計算である。

元134remainder r_iではparameter quadratic部分はh^0のQ_i(L)だけで、h^1,...,h^4の係数はparameterに一次である。monic companionでHが現れるのはtop rowへのwrapだけなので、T^j r_iのquadratic entryはtarget h^kでk≤jにだけ現れ、k=jでは正確にQ_i(L)である。producerは全134初期形状と全670×5 entryのこのpatternをQQで検査する。

四moment columns r(a),Tr,T²r,T³rのbottom h^4 rowは全て一次である。したがってbottom4 determinant delta0の次数は≤7となる。generic weighted companionでr=(L0+Q,L1,...,L4)、weight(Q)=2、weight(L_i)=weight(H)=1として計算すると、最高次数部分は正確に

    -L4*Q³

である。実sourceでQとr4が非零である証明書も保存するので、generic beta/aでdegree7は実現する。

Cramer numerator Uはbottom h^1 rowをtop rowで置き換えたdeterminantであり、なおh^4 rowを含む。従ってUもdegree≤7である。これにより以前のexact criteriaのtarget boundは

    delta0 e0 : degree≤7,
    delta0²   : degree≤14,
    delta0³   : degree≤21

へ下がる。135-generator scalar alternativeではdelta0^4 r_i(U/delta0)のdegree≤30となる（h^0係数≤2、delta0/U≤7）。norm generatorは以前どおり≤10である。

inhomogeneous ideal membershipのcofactor degreeには、これらのtarget degreeから上界を与えない。元S61のactual deltaを展開したのではなく、普遍companion determinantと実係数staircaseを検査した結果である。

保存: `code/s61_affine_normal_schur_staircase_20261103_v63.py` と `reports/S61_AFFINE_NORMAL_SCHUR_STAIRCASE_20261103_v63/`。この補足でproof packetの追加を止め、保存checkpoint後の中心課題をactual membership/unitとする。
