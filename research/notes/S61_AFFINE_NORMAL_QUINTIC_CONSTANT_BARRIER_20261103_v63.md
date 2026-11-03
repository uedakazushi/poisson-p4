# 新global quintic表示の定数cofactor障害

2026-11-03。新しいglobal monic quinticの5×670 module表示について、parameter-constant行結合では**一つも非零の一次parameter consequenceを得られない**ことを厳密に証明する。これは新しい座標と有限moduleに対する障害であり、以前の元q次数3 MacaulayやD1 left inverseを再試行したものではない。

formal coefficient ringをR'=K[n0..n41,S1,S2,H]とする。元134式のh-remainder r_iとh^j r_i (j=0..4)をmonic quinticで再び割った670個の5成分vector polynomialを全て構築する。各entryのparameter次数は≤2で、定数項はない。これは元Fがqに二次、p1,p2がhに四次、gが五次である有限division formulaから一般betaでも保たれる。

凍結QQ係数の全monic division identityを再検査したうえで、全分母が32003単元であることを確認してmod32003へ落とす。二次係数blockの形は670×5145、一次係数blockは670×225である。非零entry総数は237811である。

installed FLINTのexact finite-field rrefで二次blockの全row rank670を得た。選択670列の元のsquare matrixを作り、非零determinantも再検査・保存する。行は全670、列selectorとparameter monomial labelsも保存した。従って同じQQ minorはp-unitであり、QQおよび一般betaで二次blockのrankは正確に670となる。上界670はrow数だけによるので、char0でrankが増える問題はない。

したがってparameterに依存しない係数の非零行結合は、必ず非零の二次部分を持つ。純一次または定数のvector consequenceは一つも得られない。これはhのcofactor次数≤4を全て含み、係数はK=C(beta)まで許したgenericな障害である。

S1=z1-p1(h)、S2=z2-p2(h)、H=g(h)とすればformal presentation自体も元AとK-algebra同型である。ただしnormal42だけのbaseに戻すときS/Hは保存normal polynomial circuitsであり、新しい独立parameterとは扱わない。

この結果はparameter-dependent polynomial cofactors、より高いmodule closure、unit、一般五点尽くしを排除しない。GFでのunitを一般QQへ移した推論ではなく、最大row数の非零minorによるchar0 rank下界と自明上界の一致である。より大きなparameter次数の探索や42変数GBは今回行わない。

保存: `code/s61_affine_normal_global_quintic_rank_20261103_v63.{py,cpp}`、`reports/S61_AFFINE_NORMAL_QUINTIC_CONSTANT_BARRIER_20261103_v63/`。入力を先に保存し、native計算は60秒・2.2GB上限で行う。
