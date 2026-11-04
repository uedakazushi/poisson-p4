# 有限13次候補の基底Gに関する最終ゲート

2026-11-04。dの互除法や追加minorを計算する前に、受理済みの49項多項式G=Q L^4をKa13=QQ[a]/g13上で評価する。L=a+(9e²+3e−6)/50は4項、Qはa次数2・e次数4の9項である。e=t(a)は受理済みの全係数式を用いる。

この計算は個別のaの値を選ばない。L、Q、Gそれぞれのsourceに共通QQ分母とtの整数分母を掛けた整数多項式rawを作る。保存値n/dとの比較は、整数多項式恒等式d*raw−D*n=g13*quotientを保存して行う。元の49係数のG=Q L^4と、商環での同じ積恒等式も検証する。既知の因子以外の因数分解をしない。

受理済みの全Q/L閉成分の除外は、R=QQ[a,e,d,1/(adFBCJ)]とI=I_3(B7x3)について(I,G)=Rを与える。これはGがR全体の単元という意味ではなく、R/Iの単元という意味である。一方、norm pencilのsource恒等式、Gの単元性、既受理のKa13 source bridgeにより、Ka13[d]/H4をFULLUで局所化した環からR/Iへの全射がある。したがってG(a,t(a))がKa13で厳密に零なら、R/IではGが零かつ単元であり、R/I=0が従う。この場合はnilpotentも含めたactual Bbad schemeの空性である。

Gの値が非零の場合は、必要な場合に限ってg13とのliteral整数Bézoutを計算する。単元なら候補を保持し、真の非単元なら対応する基底成分を保持する。非零元を自動的に逆にしない。H4のd互除法、Bの別minor、Norm4、新しい点・素数・GBはこの計算に含まない。

今回のゲートの独立レビューは別途必要である。FULLU上のBのrank3が得られても、full cubic57・quartic239・generic Poisson exhaustionを自動的に結論しない。

実行結果は0.8734秒で、L・Q・Gの値はすべて非零だった。Gの値のprimitive整数分子をpとすると、保存した整数多項式S,Tと非零整数RについてS*p+T*g13=Rを全係数で検証した。したがってGはKa13全体の単元であり、H4候補coverを一部も捨てずに保持する。このゲートからactual Bbadの空性・実在はどちらも従わない。

L・Q・Gの値は各13係数で、分子・分母の最大bit長は順に1442/1408、1445/1417、1459/1444である。整数raw式はa次数24・48・144。source評価の商環演算は146回で、dは未使用のままである。Bézout多項式S・Tは13項・12項、最大20406bitで、非零定数Rは20492bitである。

再現: `.venv/bin/python code/s61_nonsquare_raw_H_B_K13_G_base_gate_20261104_v63.py --output work/G-base-gate-fresh --work work/G-base-gate-worker-fresh`。60秒・2GB・32MiBの上限で一回だけ実行する。
