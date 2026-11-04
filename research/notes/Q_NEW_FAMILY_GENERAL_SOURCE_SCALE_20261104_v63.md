# 全次数不変式に対するsource cの斉次性

2026-11-04。既存の低不変式のc冪だけから外挿せず、全16 primitive β座標と
元のSp4表現で共変性を確認する。

全source座標はβ_i=c^{e_i}b_i(u,v,X,W,B)で、
e=(-1,-1,0,0,1,0,1,0,1,2,1,2,0,1,0,1)。元Cartanの二weightを
(w1_i,w2_i)とすると2e_i−1=−w1_i−w2_iである。
実際にg_lambda=diag(lambda^-1,lambda,lambda^-1,lambda)のLie表現を
元の十行列から再構成し、その16固有値が2e_i−1と一致することを検証した。
従って全source上で、単なる有限点・特定wordでなく

    beta(u,v,lambda² c,X,W,B)=lambda rho(g_lambda)beta(u,v,c,X,W,B)

が成り立つ。Sp4不変な斉次degree d多項式なら、cについてc^{d/2}の単一冪である。
中心−Iにより奇数dは零。これはd=14,18,26等にも使えるが、実際の新しい不変式
展開は行わない。既存section vはcに依存しないので同じ恒等式が保たれる。

xi=40c/(Bu)を用いる既知involution固定source体では、14次不変式は
xi^7 h(shape)となる。定義次数下界26により、非零の普遍14次不変式はNEW上でも
非零だからhも非零である。

**条件付きの次手。** spectral比とZN=N12/F4³が全shape体を生成すると別に証明
できれば、F4=xi² f(shape)とI14/F4³=xi*h(shape)/f(shape)³からxiを回収できる。
このときactual target不変体と既知固定source体のindexは1になる。このpacketは
ZNのprimitivity、shapeの回収、target index1自体を証明したものではない。

再現は同名Pythonへ新しいoutput先を渡す。一回60秒・2GB上限。元の全16有理式の
c除去と、標準4×4から元16×16へのLie表現を全係数で照合し、証拠とhashを保存する。
