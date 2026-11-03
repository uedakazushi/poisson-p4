# Q_b actual D cubic gateの独立物理監査

2026-11-03。対象は `AMBIENT_CB_ACTUAL_D_CUBIC_SYZYGY_20261103_v63`。
full D theoremそのものは別applicationとして監査する。
元packetには最終b/c/cofactor arrayとinverse/Jacobianがあるが、raw source offset jetは保存されていなかった。
そのため「保存されたraw source jetを照合した」とは主張せず、独立に必要な八曲線を再構成して
本REVIEWには全offset jetと全physical outputを保存する。
旧82 callback、22秒worker、98秒SymPy workerは再実行しない。

## 別のphysical D3 evaluator

全150係数の外積indexと二次monomialの逆座標微分から、各Eijの疎Lie作用素を独自に作る。
author helperのevaluate/evaluate_lieは使用しない。
CA20とCBのf/H/ratio sourceを固定し、CBでは全degree0–3のmu*f convolutionを保つ。
flag offset `A(t)=A1*t+A2*t²+A3*t³` に対し、三次係数を

    F3+L1F2+(L2+L1²/2)F1
       +(L3+(L1L2+L2L1)/2+L1³/6)F0

として計算し、最後に実固定g0の逆pushbackを施す。
行列squareと150係数Lie作用のsquareを混同せず、L1L2の順序も保つ。
全三sourceの150×29 original Jacobianを保存値と照合する。

実58/87 constraint関数を全物理係数から組み直し、全square matrixとinverseをQQ(b)で検証する。
一次解から二次・三次のknown source residualを独立式で求め、未知degree-d offsetが
元Jacobianでlinearに入ることを使って補正する。
mode58/87それぞれ四つのslopeで八曲線を作り、全constraintのdegree1–3を確認する。
normal三次は四monomialしかないため、slope 0,1,−1,2の明示Vandermonde係数式で完全に復元できる。
bは最初から最後までQQ(b)のsymbolic parameterであり、b標本からの零外挿ではない。
0次冪は明示的に1とし、0^0のadapter例外も避ける。

## mixed cofactorと全150/121 coverage

七つの残りlinear conductor方向とx0/rhoの組に必要な十四mixed Hessianは、
源のmu*f mixed項および対称化Lie二次式から直接計算する。
各々を58 constraint inverseでprojectするので、三本のsource callbackを二十一回繰り返す必要はない。
これで全150のH0、H8、Hiを独立に復元する。

common base29を除く全121 physical座標を保持する。
実f row18のc0=1、他全149行のconstant c=0、quadratic generator row63のc_quad(0)=0を明示確認する。
linear generatorのcofactorはdegree2、quadratic generatorはdegree1まで必要であり、
quadraticの未知degree2 cofactorはc_quadの正次数と掛かってdegree3以上となるので目的のquadratic defectへ寄与しない。
他の七linearのc_iも正次数で、cofactorの未知degree2は同じく目的次数に寄与しない。
CA section上のそれらのb_iはdegree2以上なので、全b-syzygyの三次確認にも計算したdegree範囲で足りる。

全150のb-syzygyを三次まで独立に検査する。
全c-defectのconstant/linearが零、全quadraticがactual Z4の二quadric spanへ入ることをQQ(b)の
完全係数等式で確認する。従って唯一残り得るcubic socle classは全150行について恒等零。
これは数個のambient rowだけやb=1だけの検査ではない。

本gateが示すのは全ambient redundant odd vectorを含むcubic障害の零である。
高次normal冪の既知core idealへの包含と相対free rank1 moduleへの帰着、dominant G×bによる適用は
別full D proofで明示する必要がある。このgate単独でnilmodule/Ass/global成分を宣言しない。

保存先：`reports/AMBIENT_CB_ACTUAL_D_CUBIC_SYZYGY_20261103_v63_REVIEW/`。
