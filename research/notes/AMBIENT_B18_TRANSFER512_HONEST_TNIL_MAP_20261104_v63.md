# 全512 transfer cycle の honest Tnil128・minimal5 への実写像

2026-11-04。元の実44-cycle全512を ell=0 に制限し、canonical m24×44 を掛ける。各原 column の全24 formal coefficient functionに三 sector の ordered Hadamard・Euler homotopyを適用して b84 を得る。d2b=mz なので128-vector(z,b)は Phi_nで zへ戻り、実 W→N/ellN→Tnil を与える。先に係数をfactorしたり有限jetをhomotopyへ入れたりせず、全512個を構成してからS-linearに延長する。全44、128、5の演算子を有限 typed DAG に保存する。

## 1. 二つの128表示を混同しない

元 honest128 は section n、minimal演算子の base128 は E-supported section n′を使う。J_n(s,b)=(Phi_n(s,b),b) と逆側J_n′を明示する。両者は actual(X,h)、section、d3関係を保ち、合成とidentityの差は対応section-boundaryである。512個の選んだliftはPhi_n(z,b)=Phi_n′(z,b)=zを満たすので、J_nが実際にその128-vectorを固定する。従って既存minimal5×128 projectionを適用できる。二つのpresentationを単に同一視しない。

## 2. cacheだけによるorigin検査のsource gate

全512 actual W44のoriginは零である。a36の全項にβ_iが入り、odd cの各項にもb_j,b_Aj,δ_jのいずれかが入り、これらは同じoriginで零。実C4 cache、actual A0、raw physical F_j、C0 source第一係数を使い全25,088 typed source nodeのoriginを評価する。

最小projectionの残るwedgeはsector0の(0,7),(4,7),(5,7)。m0row7=(5/2)e20であり、W0=0なのでm1寄与は零。従って必要なのはW20のz0,z4,z5一次係数だけである。actual sourceは

    W20=even_C[β_i A_j,5 (p_B+τA p_AB)+indicator β_i b_j e_gamma]。

全p=n_alpha e_gammaは原点でp_B=0。nil n1のoriginはτAτB成分のみ、n2は零である。C/m_R CではτA²=τB²=0なので、任意のactual A0についてeven(A0 τA)=0。β0=0のため他のproduct-rule項は零となり、W20一次係数は全て零である。未知source/gauge derivativeを零に置かず、零係数に掛かることを確認する。

## 3. 結果と停止範囲

全512のminimal5 origin columnは零、rank零となる。したがって imW⊂m_S Tnil、μ(Tnil/imW)=5、長さは少なくとも五である。強い十分条件「Tnil/imWの長さ二」は失敗する。

これはWの像そのものが零、N³非零、またはP=psi(W)が偽という結論ではない。正しいP被覆用のK/ellK注入・長さ九moduleと正規化fiberの不可視二〜四方向を保持する。全scalar512 homotopy係数や新source・jetを展開しない。
