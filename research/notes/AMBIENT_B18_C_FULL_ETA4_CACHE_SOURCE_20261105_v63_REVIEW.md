# 四つの新しい eta 係数の既存 cache 読戻し：SOURCE・実装監査

独立報告 `f18e258a...` は主 SOURCE `5052f2af...`、固定 worker `c6e00c67...`、
custody `eb37e44b...` を実値前に監査した。親の条件付き一回実行条件はこの PASS で満たした。
本監査自身は要求された正次数 root や eta を計算していない。

四行19、28、29、72の八 root は originalID を持つ7945演算の最小 ancestor closure と一致する。
同じ実105 unknown・105 equation の既知単元 Jacobian 証明書を保ち、
r=0、s=ell/5 の環準同型で提供済み次数三係数を移す。
105入力の既知315係数と八 root の既知24係数を対応付けた。
これは旧 alpha=delta_r ell² とは別の、新しい eta=delta_ell³ の source である。

45個の分数定数は分子・分母を両方保持し、32個の非整数を確認した。
正の絶対上界だけを用い、逆数は独立の閉じた三次公式で評価し、
積は逆順に加算して全7945演算の上界 hash を再構成した。
二つの逆数検算積と delta 差の全 partial も含め、保存上界35bit、演算前71bitである。
要求値の相殺を上界に利用していない。

worker の分数 parser、各 add/mul 前の2048bit guard、逆数 recurrence、
両側の積検算、全 lower overlap、last-use release を文字通り確認した。
live node は184、独立 SOURCE 全体上界116069、将来値 worker は peer5000枠込み55985で、
静的75000・hard100000・全体150000を下回る。RAW は別に管理し32MiB以下である。

結果の数学的範囲は、完全155 source の四 eta による全 C の HF3 五又は六の決定だけ。
この worker 自体は長い DVR 指数、元 H0、nil 三乗を決めない。
後続の条件付き純粋定理を適用する場合は、その別の証明と値検算を必要とする。

再現は `code/ambient_B18_C_full_eta4_cache_source_review_20261105_v63.py`。
