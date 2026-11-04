# κ・tail2 actual11pivot の独立 source／演算子監査

2026-11-04。対象は projected prePOWER15×22 とその全 tail である。新しい
parameter coefficient、8minor 展開、source jet、rank search は計算しない。

先行する独立 QQ peer は全238 source 値、四 signed cofactor、B を掛けた追加列と
四積を検証し、T の7208項非零を証明した。anchor は文字通り M[:,22] であり、
余分な a はない。本監査はこの peer、元六pure列の rank≤1、実70 heavy係数の
known-unit frame、Saux の既存 source と非零証明を固定する。

四 bottom source の κ は [Hbar0,Hbar1,M22] の cofactor である。元六pure列の
rank≤1 と非零anchorから、κ は既知pure列も多項式恒等式として消す。anchor を
各点で反転する議論ではない。最後の行を κ 結合に変えると determinant は κ3 倍、
唯一残る列は index5 の T である。符号 (7+5) は正。

七行の補小行列について、三つの bottom row の係数 c1,c3,c4 を抜く。
先頭四行の heavy二列が零なので、その二列に沿う一般化Laplace展開は

    det7/(c1*c3*c4) = det5(P4, δ2*Gu/c4 − μ0*G1/c1 − μ1*G3/c3).

右辺は保存 source の Bhat と bordered Saux の有限式そのものである。
独立な小さい自由可換多項式 checker で2×2 adjugate符号、三行Laplace式、
全768 source式の cut版 bordered式を照合する。heavy行・列の移動は8回・10回で
ともに偶数である。この証明は多項式のまま行え、δ2 や p の反転を要しない。
Saux の7022有理数係数全体を新しく再展開したとは主張しない。

従って κ3 detA8 = cprod Saux T、cprod=−5120000000000000/45927。
D(Saux*T)では積が単元なので detA8 の逆は κ3/(cprod Saux T) である。
κ3 を別に反転しない。三つの追加B source は current15 [7,13,10]、heavy列は
元32 [28,30,29]。前八行のこの三列は literal zero で、右下blockは
diag(rho6,rho12,rho9) N3。Adj3(N3)/Delta407 による逆の順序も照合する。

全84 typed node、外部hash、source／target permutation、列2のB倍、heavy二列の
L=aFC² と L gamma による除算を検証した。自由非可換block代数で K11 の両側逆、
Schur左右演算とそれぞれの逆を検証する。旧operatorから参照するのは実入力、
B3およびcurrent15→old25 source prefixのみで、旧H5逆を持ち込まない。

受理する領域は **FULLU∩D(Delta407*Saux*T)**。同cokernelの残余は4×11、
sourceはcurrent15 [9,11,12,14]、targetは元32 [1,3,4,6,7,8,9,14,17,18,19]。
tail2をpivotとして消費し、他七tailと全補集合を保持する。p、R434、h22、Delta5
は反転しない。raw H7 の eta/HZ/HW/Cdelta lift は別gateであり、残余rank4、
全FULLU57、全quartic239を結論しない。

再現: `python code/s61_nonsquare_prePOWER_kappa11_schur_review_20261104_v63.py`。
