# 個別 native pivot 二次九十二係数 worker の準備

この worker は immutable SOURCE ab2e の元23列と current Gamma56 を保ち、元 N_B11(rr)・N_C16(ell/r ell) と canonical C68 の G8・monic だけを使う。N転置は十二源・十二物理行だけで、四trace blockを転置しない。元 native X の下八行は零であり、Gamma自体を48行と取り違えない。

Gamma Res4(N転置) の五plane係数を内部で求め、Res4(S) の各23列と ordinary coefficient の積を作る。rrは三項、r ellは四項で、最後の L二次 S0 を全て保持する。free16列はこの最後の項だけが寄与し、weighted-kernel の S0零を個別列に用いない。tauのcollisionに伴う N0・nu二次とintrinsic Sのmonic carryを両方残す。

出力は二channel・23列・二pairの92係数と、検証用の322個の積則内訳だけである。SOURCE280 committed・841 temporary・全保持21034を維持し、失敗時のcursor・GammaL prefix・intrinsic S・完了値を保存する。新しいphysical/gauge解、honest quotient、N三次、617 upper、HF・全長の応答は求めない。準備時点では値未計算であり、ROOTの明示GOと独立SOURCE peerを必要とする。
