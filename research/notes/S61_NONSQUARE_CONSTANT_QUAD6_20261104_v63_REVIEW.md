# Constant quad6の独立proof audit

2026-11-04。全raw aC開で六quadの定数階数が6、actual aFBC chartでweight2商dimension23という新しい結論を受理する。22 quadratic syzygiesの独立性をd=0/J=0へ延ばす主張は受理対象でない。

独自Fraction dictionary arithmeticで、保存されたactual ell/Lから三つの必要な係数恒等式を再構成した。ell_x1=−16CJ/243、3Delta_x2−(7/2)Delta_x4=−80ad²/81、d=J=0すなわちd=0,e=−2/3でell_x0=3200a/243である。authorのdecoderや実装をimportしていない。

aC開でJ!=0、J=0かつd!=0、J=d=0の三場合は互いに排反かつ網羅的で、各場合に対応する係数は非零。従ってellは零多項式にならない。d=0かつJ!=0は第一場合に含まれる。

既受理のZW20×35 rank19から四quadが定数独立となる論証も正しい。非零定数関係cがあればcを五つのlight変数で掛けたsource vectorsは五つ独立（v↦c⊗vがinjective）なので、kernel dimension1に反する。

元六quadから(q0,q1,q3,q5,ellZ,ellW)へのdet−80、分母aだけのbasis changeを使う。六式の定数関係をZ=W=0へ制限すれば最初四係数が零。残りell(bZ+cW)=0はpolynomial ringのdomain性・ell非零・Z/W独立性からb=c=0。これでraw全aC上rank6である。

actual weighted chart aFBCでは元式とのunit row/frame変更を戻せ、weight2の全29単項式に対する六独立関係から商dimension23となる。高weightの式はweight2へ寄与しない。full cubic57、quad146、fullquartic、ambient除去には延長しない。

新しいrank・field・point・GBは計算せず、全caseと既知unitの論証および少数のactual coefficient identitiesだけを監査した。
