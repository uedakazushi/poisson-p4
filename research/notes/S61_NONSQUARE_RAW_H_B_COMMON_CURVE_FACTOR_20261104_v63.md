# 残存共通必要曲線Gのexact factorization

凍結済み49項Gだけを因数分解する。a/e次数6/12、primitive整数係数39bitで、
入力のmonomial boxは91以下。native FLINTのQQ multivariate factorと、別の
SymPy QQ factor_listの全monic因子・重複度を照合する。全係数の整数分母払い
再構成も確認する。新たな点・素数・GB・終結式は使わない。

出力因子は共通根曲線の必要成分であり、実際のB-bad点の存在を意味しない。
次に各成分でH,g,kのpseudo-remainder又はsubresultantを調べるなら、分母が
成分上で恒等的に0でないことを確認し、消失する有限残余は別に保持する。
FULLUに入らないJ/B/dの境界で全候補が説明されるかはまだ未検査である。
