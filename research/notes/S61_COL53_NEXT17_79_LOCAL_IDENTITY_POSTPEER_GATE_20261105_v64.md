# 保存セルの局所恒等式による独立 postpeer の準備

これは primary `4b1ed027…` の実結果を検算するための**未実行の独立 signed postpeer**の SOURCE ゲートである。prepared report は `reports/S61_COL53_NEXT17_79_STREAMED_VALUE_20261105_v64_LOCAL_IDENTITY_POSTPEER_GATE/report.json`、SHA-256 `9eb3ed13ecab356c8001366e636d769e27def438a299aa8bd7500367e9034e95`。worker は `code/s61_col53_next17_79_local_identity_postpeer_20261105_v64.py`、SHA-256 `492200c950fdcdf6e2d6894f3f5dcc6c41cb421934cce7c17803c3f9b06298c5`。準備段階では符号付き node 演算を行っていない。実行には、この gate と code を直接 pin する ROOT の新しい GO receipt が必要である。

検算法は生産者の全 prefix を再評価する方法ではない。manifest を一行ずつ読んで offset/bytes/count/hashだけの thin index に縮め、各 node の保存済み current・保存済み children と真の元 leaf に対して COPY/ADD/MUL の**局所係数恒等式**を別実装・逆順の項演算で検査する。1884の全恒等式が一致すれば、topological inductionにより root と target の証明になる。各 temporary canonical cell と immutable stream segment の二世代をバイト単位で比較する。生産者の係数 evaluator は import しない。

独立 bounded parser は stored primitive content と固定 source q_i を区別し、G=stored/q_i を正整数として c_i=d_iG を復元する。空の cell は9バイト flag0を要求し、q_iによる除算をしない。primitive整数、分母、項数、canonical順序、gcd、5座標、M_i、SOURCE cardを検査してから係数を使う。元の leaf も literal unit/monomial transform、既存 binary payload/hash、元の full5座標から結び付ける。

各保存 child が局所恒等式で検査済みであることを使う。固定 q_i の copy/add/product weight は整数で、元 SOURCE の正値 partial bound≤M_i を準備段階で全1884について再確認した。整数係数824ビット、復元・一時演算1921ビット、hard2048を保持する。恒等式の差も高々2M_iであり、新しい高位係数や parameter point を求めない。

数値保持は full manifest を避け、thin index と一つの current、二つ以下の child、独立 expected mapだけにする。各 map 項は全5座標と係数・slotを7数値として数える。全1884の検査記録とその serializationを含む metadataの保守上界は107,607<META118,000。元 leaf は一つの literal term表現だけを作り、全-key mapからliteral表への二重コピーをしない。元 SOURCE の decoded14104上界から最大2349項であり、入力 transient21000を割り当てた。

SOURCEの card41を持つ入力で、current oracleとexpectedを同時保持した通常上界は143,670、partial保存を含む停止時上界は**143,793<150,000**。他の codec/decode/root-copy phaseも固定SOURCE全1884から計数した。二つの canonical世代とmanifestを含む元生成物23,968,000、既存leaf入力、独立出力reserve等を含むRAW上界は27,389,132<32MiB。AS2GB、55秒signal、60秒CPUの範囲を維持する。

準備時に実際の保存係数を演算したのはなく、node fileの hashとサイズの tree digestだけを固定した。無関係な既存 synthetic codec例で bounded current parserと単一literal parserを確認し、後者は別のencoderによる四つの厳密roundtripも確認した。draftの構文検査で見つかった説明用dictのcolon不足は、prepare前に snapshot/historyを保存して修正した。生産者・primary・以前のgateは変更していない。

成功しても結論は選択した H/L固定・relative s≤2 の129候補射影とtargetに限る。whole minor の零、rank、次のpivot、matrix/RHS更新は主張しない。STOP時には検査済み局所恒等式のprefixと現在のexpected部分を保存し、自動再試行やcap変更を行わない。
