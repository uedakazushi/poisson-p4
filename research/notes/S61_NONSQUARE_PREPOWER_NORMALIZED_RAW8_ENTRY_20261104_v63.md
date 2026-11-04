# 固定normalized raw8 determinant：二重整数融合の一回計算

2026-11-04。入力は凍結済みNESTED_RAW8_SUPPORT_PREFLIGHTの同じ一候補のみ。
source rows current15 `[0,1,3,6,2,4,5,8]`、元32 columns
`[22,11,13,25,15,19,27,31]`。FULLU既知単元L=aFC²とL*gammaで二heavy列を割り、
既存証明のa row冪だけを除いたdeterminantを計算する。

全source原子・1898通常演算・90inner/15outerの理論QQ分母を固定する。bufferは全degree/box
格子点を持つ58259 outer、一つだけのbottom≤14484で、各畳み込みをfmpz_addmulで直接
蓄積する。実contentから得る整数scaleも全eventで検査する。新polynomial productや
accumulator全copyは作らない。通常source演算のscratch、preloaded atoms、integer分母と
work reserveを含むガード150000を維持する。

各event後にprogressを保存し、50秒nativeガードでは完成済みouter項数・inner項数・途中
畳み込みの左/右index、exact partial outer/bottom polynomialを保存する。この場合は
未完成の行列式や非零minorとは呼ばない。別candidate、prime、factor、再試行へ進まない。

完成した場合のnormalized polynomialに、`L²*gamma*a^(sum rowa)`を掛けたものが元の
actual8×8 minorである。この既知単元を未展開のまま保持する。非零であれば、既存B3を
加える11pivotのdetはそのminor×Delta407×既知rho単元となり、指定chartでsame-cokernel
4×11へ進む。値が零又はguard停止なら、そのままinconclusiveとして保存する。
いずれでも全FULLU cubic57や残余rank4は結論しない。

sourceは全Q9tailを含む元15×22からの正しい部分minorであり、将来の4×11にも全8tailを
保持する。p、R、Saux、h22の逆をこの実計算へ導入していない。
