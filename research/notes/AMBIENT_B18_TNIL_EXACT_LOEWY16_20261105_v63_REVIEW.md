# 実 filtered16 basis の独立証明

五つの新minimal generatorと一次の五単項式を実sourceに照合し、二次の三列 `i(e3*r),i(e4*r²),free42*r²` の固定行列式を確認した。最後の列は文字どおりの grQ2 双対に検出され、前二列とは独立である。

`u=e3*r,w=e4*z5,q2=e4*r²` とすると、未知alphaに対して `u+c*w-alpha*q2` は旧moduleのr核に入る。その新gr2像は `(1+c*a)i(u)+(c*b-alpha)i(q2)` であり、第一係数が既証の非零有理数である。r:gr2→gr3の核は一列なので、これと `i(q2),free42*r²` の行列式はalphaによらず非零となる。従って `i(g),free42*r³` が三次の基底である。

四次は実source式 `d=i(qC-a*e3*r²-b*g)` が生成する。qC係数が非零で、旧alphaの値は不要である。これにより全16 liftはfiltered basisとなる。Hilbert関数は(5,5,3,2,1)、冪の長さは16,11,6,3,1,0であり、m²の k[[rA]] 分割は(3,2,1)である。

完全917 operatorのrank1とright106の零からconstant811のrank1も従うが、特定の非零列や全作用行列は計算していない。rAはDVR ellとは別であり、このbasisは元nil idealの生成元や積を与えない。
