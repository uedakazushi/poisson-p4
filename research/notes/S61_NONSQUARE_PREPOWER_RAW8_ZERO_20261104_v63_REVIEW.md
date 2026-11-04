# 固定 raw8 の零性：独立 QQ 監査

対象はrows current15[0,1,3,6,2,4,5,8]、columns original32[22,11,13,25,15,19,27,31]という一つの候補だけ。新しい点・prime・rank探索・他列を計算しない。

既受理のsource/lifetime gateから下段pure4×6に必要な190 node closureを再構成する。保存された1898worker inputと190抽出inputの全operation/参照を比べるだけでなく、46 literal atomの全QQ係数を元source JSONから独立に再読し、aの分母/shiftとQQunitを一致させる。四行のrow-a除去は全て0であり、得られる24値は元actualMのpure六列そのものである。

独立native checkerはproducer headerを使わず、異なるDEGREVLEX monomial orderで190個の固定source演算をQQ上で再評価する。24個全てを保存値と比較し、その後90個の2×2 minorを二つのQQ多項式積として直接照合する。大きいnested raw8 workerや出力accumulatorを再実行しない。

90個が全て零なので、この下段pure4×6は任意の幾何点でrank≤1。下段の残りheavy二列はrank≤2で、下段4×8のrankは≤3。上段は四行でheavy二列がliteral零であるため、全8×8のrankは≤7。従って選択raw8 determinantは恒等的に零。既知単元L²gamma*a²の正規化もこの零性を変えない。

この候補による11pivotのchartは得られない。四行の他のpure列、全17pure列のrank、他候補、全cubic57やquartic239について否定的結論は出さない。未知の因子を単元としたり、零minorから実欠損locusを推測しない。
