# delta63 の一つの選択的四次 functional：origin adjoint だけを固定

残る delta63_r²ell² に対し、同じ原105方程式・未知数順序と固定単元 Jacobian を使う。
原ノード20144−40922を凍結済み原点で reverse AD し、
rho=(grad_Phi delta63) J105^{-1} を得る。保存済み J と逆行列の両積、および
rho*J=grad の全105成分を厳密に検査する。新しい逆行列や IFT を解かない。

次の値計算では B9=QQ[r,ell]/(r³,ell³) を使えば四次単項式は r²ell² だけ。
提供済み全平面三次 Phi を入力し、四次 slot を零にするのは lower-term と
方程式 RHS の組立てだけであり、真の implicit Phi4 が零という仮定ではない。

真の値は g_lower4−rho*RHS4 となるが、本 packet はその値や RHS4 を計算しない。
rho が読む方程式と二つの g root の原 DAG closure、last-use 順だけを保存する。
全 node の正の高さ、coherent post作用の予算、prepared code/custody、独立 peer と
親の新たな ONE GO を揃えるまでは四次応答を問い合わせない。

g は even difference delta63 であり、odd f63 や gA と名前を混同しない。
元 H0、K torsion の指数、nil 三乗への結論はない。
