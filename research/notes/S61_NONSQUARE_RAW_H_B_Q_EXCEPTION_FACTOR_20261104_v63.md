# 安全なQ例外37次式の厳密因数分解と成分inventory

分母例外を全て保持して得たSF37だけをQQで因数分解する。38項、primitive整数
係数198bitなので小さい入力。native FLINTと別実装QQ factor_listの全因子と指数を
照合し、元の37次式を全係数で再構成する。点・素数・GB・終結式は使わない。

各因子のe次数mに対し、Qがaモニック2次、Hがdモニック4次なので、a局所化前後の
有限QQ代数は次元高々8m。matrix N×22Nのentries上限も記録する。
最小次数・項数・native factorindexの順に一成分だけを次のmembership候補に選ぶ。
他の全成分は未検査のままcoverに残す。因数分解自体では空性を主張しない。
