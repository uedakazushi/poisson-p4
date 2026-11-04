# 有限例外のsquarefree化と小さな有限代数入力

L/Qの保存済み例外polynomial（積次数65/361）だけを使う。QQ単元でモニック化
して積Pを作り、微分とのextended GCD G=sP+tP'、S=P/Gを厳密に検証する。
SとS'にもBezout=1を保存する。さらにPからSを割る因子だけを反復除去して最後に
1になるliteral quotient列を保存し、PとSの零点集合が同じことを全係数で確認する。
点・新素数・GB・終結式や有限fiberのunit検査は行わない。

squarefree次数nL/nQにより、L側の候補は高々4nL点、Q側は高々8nQ点となる。
a局所化前のbaseはL側QQ[e]/S_L、Q側QQ[e,a]/(S_Q,Q)で有限次元。
aの逆元を加えた後もQQ上の有限代数のdirect factorであり有限、さらにモニックH
でd方向に高々4倍。これはmodpへのproper integral modelを自動的に与えるもの
ではない。古い有限残余の高々1600点を加えた幾何学的上限だけを更新する。

oldh1のe16最高係数もliteralに確認する。これはa多項式ではなく非零QQ定数。
従ってh1をeモニック化でき、Hと合わせQQ[a,a⁻¹]上の有限free64モデルになる。
k1=e25の最高係数もQQ定数として保存する。ここではk1によるfiber検査や
新たなnested resultantは行わない。B-bad空性、globalB3、p=0/full57/239は未解決。
