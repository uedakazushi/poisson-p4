# Higher Bianchi の一変数除算: generic parameter-field gate

2026-11-04。既存 `GENERIC_RANKS_20261102` の同じ点 (a,e,d)=(1,2,3) だけを用いる。新しい β 点・field・source 計算は行わない。

凍結済み weight5/6/7/8 行群の20/22/18/11本を、今度は QQ(a,e,d) 上の未知係数で結合する。factor候補 v の weight は w−4 とする。全135 cofactorの v で割れない係数を `(generator, local multiplier variable)` ごとに集めた **unsplit** matrix を使う。a,e,d の係数を別の constraint に分離しない。19/23だけは実 chart の Tnormal/U 両行が literal zero であることを再確認して代入する。

各 matrix を同じ既存点で有理数として評価し、最大22次の square minor とその **QQ の正確な行列式**を保存する。この値が非零なら対応する polynomial minor は恒等的に非零であり、QQ(a,e,d) 上でも columns は独立である。従って同じ weight 群の非零 parameter-rational 結合を v で割ることはできない。点で rank が落ちた場合は generic rank について結論しない。

43 systems 全部について、選択した全行・全列・有理数 matrix・determinant を保存する。F5/F6を含む135列を保つ。これは source-level の指定した exact-weight4 一変数除算 class のみを扱う。別の weight の変数によるより低い weight への除算、local variables の線型結合による因子、一般の gcd、saturation、normal graph 上だけで可能な除算を排除しない。全点での injectivity、新 pivot、全quartic rank の主張もしない。

一回だけ、上限30秒・500MB。元 operator の再生成・GB・新しい field/point は使わない。
