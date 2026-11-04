# P100上のfull cubic一般rank57：独立監査

著者のrank workerは実行せず、P80で独立検証済みの一般QQ polynomial checkerをそのまま用いる。新しいreview driverはP100の有理写像を独立にbindし、元checkpointの全408個の有理成分を再構築する。C++ checkerの算術はP80 reviewの凍結版と同じで、行列係数の期待最大bit数だけP100の153にする。初回の169というmetadata設定による停止は、全408成分一致後・determinant前であり、元sourceとlogを別failure packetに保持した。

P100ではa=−10(zeta+20)/(9zeta二乗)、e=(zeta+100)/(9zeta)、d=−3zeta(zeta+20)/100。これらを受理済みTheta有理座標の式と比較して確認し、元のB多項式のpullbackも全係数で検証する。全分母は非零QQ定数、a、Bのいずれかであり、P100のpermitted open上の単位である。元17対24行列と保存されたrow clearingを全て検査する。

固定された17列のdeterminantだけを独立な固定順Bareissで再計算し、全exact divisionを掛け戻す。非零なら残るrank17、既存八unitを加え25、既受理の元source同値でfull cubic rank57がこの曲線の一般点で成立する。

これは全P100でのrank57ではない。保存minorの零点に有限例外が残り得る。他曲線・full quartic239・全分類は結論しない。新しいbeta試験点、素数、rank search、POWER cofactor展開やGBは使わない。

再現は `code/s61_nonsquare_Theta_P100_cubic_generic_review_20261104_v63.py --output <新出力先>`。既存Cloud activateでFLINTのcompile/runtimeを再利用する。
