# 旧有限残余のeモニック16次ノルム入力

h1のe16最高係数とk1のe25最高係数はいずれも非零QQ定数。QQ単元でモニック化し、
QQ[a][e]/(h1)の基底1,e,...,e15でk1を掛ける16×16行列を作る。全16列の
k1*e^j=h1*qj+ΣMij*e^iを独立整数分母払い演算で確認する。
未知のa-leading polynomialは逆にしておらず、その消失場合も作らない。

行列entriesはaだけの多項式。整数分母払い後の係数上限と、実際のentry次数行列に
対する最大重みassignmentから行列式次数の上限を保存する。このassignmentは
整数の組合せ上限であり、行列式そのものを計算したものではない。

既存h1/k1の互いに素性とGauss補題により、提案ノルムは非零。これを将来計算すれば、
有限残余の必要a値のunivariate polynomialが得られる。ただしこのpacketでは
16次determinant/adjugate、nested resultant、GB、点、新素数は計算しない。
QQ[a]上の有限free16モデル（Hを加えるとa非零でfree64）を入力として保存する。
