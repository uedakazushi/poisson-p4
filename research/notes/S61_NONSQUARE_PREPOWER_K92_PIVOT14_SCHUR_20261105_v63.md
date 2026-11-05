# K92 の第十四 pivot と source 忠実な一対八

2026-11-05。これは新しい三変数係数展開を伴わない有限作用素である。
軸上 Bareiss primary の非零性を使うが、その独立確認は別 gate として保持する。

旧 PIVOT13 の二対九行列を Y、行を current15 の [9,14]、列を original32 の
[1,3,4,6,7,8,9,17,18] とする。Q9 の他の tail、消費済み tail2 の効果と
full57/canonical63/pure84/raw7 の source は、旧作用素から全て引き継ぐ。

P_K92 は保存された正規化済み十一角行列 M の行列式そのもので、展開しない。
M の列は [22,11,13,25,15,19,14,18,28,30,29]。最後の三列の既知単元 a F C²
による正規化、G18/G39 の a 正規化、全 H3 項と符号は source recipe に固定する。
旧 fixed det11 = a² P_K92 である。従って Y の (0,8) 成分 s は

    s = c8 g P_K92 / (Delta407 N7),  c8 = -88000000/13851.

ここで g は raw_compact14x5[6][1] で、G20 の anchor ではない。
既存 N6 chart の積単元恒等式によって g は単元である。
N7 は既存の非零有理数倍 a九乗 N7bar なので、s の逆には

    Delta407 N7 / (c8 g P_K92)

を使う。新たな g、D5、p、R434、h22、generic-cleanup deltaP/deltaQ/N0 の
開条件は加えない。新しい開条件は P_K92 のみである。

Y の列を [8,0,1,2,3,4,5,6,7] に並べ [s b; c d] とする。
左作用素 [s逆 0; -c s逆 1]、右作用素 [1 -s逆 b; 0 I8] とそれぞれの逆を
明示し、残余 d-c s逆 b を得る。free Laurent ring における左右逆と
block-diagonal 恒等式を全係数で確認する。パラメータの新係数は評価しない。

fullsource15、fulltarget22 とその逆を旧作用素に合成する。残る source は
current15 [14]、列は original32 [1,3,4,6,7,8,9,17]。
残る一行の元25/full57/canonical63/pure84/raw7 source は全て
[-c s逆,1] を旧二行に掛けたものなので、raw7 の零化と実 cubic への持上げを保つ。

旧 Omega = D(Delta407 Saux T N6 N7bar) 上で P_K92 を追加した chart に限る。
この作用素は prePOWER15x22 と同じ cokernel を持つ。残余一行が全列を張るなら
既存 FULL57_SOURCE_LIFT により実 cubic57 が従う。この十分条件を、実 cubic の
不良 Fitting ideal と同一視しない。残余一行の rank、全補集合、全 FULLU57/239、
一般五点尽くしは未証明である。

軸 e=d=0 は FULLU 外であるが、全 generic 指数が非負な多項式 M への係数制限である。
primary の [a十二乗] detM(a,0,0) は非零。従って P_K92 は非零多項式になる。
この利用の独立確認が終わるまでは第十四 pivot の非空性を条件付きと明記する。
