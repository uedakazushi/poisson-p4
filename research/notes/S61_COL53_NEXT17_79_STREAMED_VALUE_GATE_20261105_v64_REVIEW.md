# 固定17/79射影のストリーム VALUE コードに対する独立 SOURCE 査読

worker `a946d904…`、gate `c5823b8d…`、custody `e75cf2b3…`、phase companion `2577d4a1…` を査読した。独立報告 `reports/S61_COL53_NEXT17_79_STREAMED_VALUE_GATE_20261105_v64_REVIEW/report.json` の SHA-256 は `012f7d578a5240351e40451479f4caeb22997b85247a9e2956de00f572cab27a`。**SOURCE・コード・custody の監査は PASS**。固定17/79の符号付き VALUE は実行せず、ROOT の別の実行許可を待つ。

既に ROOT が確認した2998-node元式と射影 SOURCE をハッシュで結び、needed1884ノードについて独立に全整数復元係数と正値 partial bound≤M_iを検査した。copy/product の5座標 shift、各子の needed relative-s degree、和で残す H/L-tied 枝を確認した。root の正値支持は129個の異なる**全5座標**で、target はその支持に含まれる。実係数の相殺が起きても source の lead や q_i を更新しない。FULL F11 の元 root1202を411項の H-faceで置き換えていない。

各 node は固定 q_i と整数係数 c_i による q_iΣc_ix^i として演算する。copy/product/add の係数比が整数であり、全中間部分和・積が M_i 以内であることを1884ノードすべてで確認した。全5座標を保持するので、同じ積の座標へ寄与する項数は少なくとも source 支持の min-cardinality 上界で制御できる。出力を primitive cell にする際だけ G=gcd|c_i|を取り、stored content=q_iG、primitive body=d_i=c_i/Gとする。子を読む際には **G=stored content/q_i を正整数として復元し、c_i=d_iGへ戻す**。stored contentを次の演算の q_i として使わない。負の q_i でも G は正で、stored contentの符号は q_i と一致する。

空の射影では mathematical G=0、codec sentinel G=1を区別する。9バイトの `Q5CELL1`/flag0だけを保存し、q_i=0の場合も q_i による除算を行わない。非零 q_i で実射影が零になった場合も、それを元関数全体の零とは呼ばない。source lead を actual support に移動することもない。

実 codec writer/validatorに対し、別実装の encoder/parserを用いる**無関係な合成例だけ**を試した。負の q、混合符号、primitiveでない固定-q係数、単項、負の座標、空の q=0/非零q、および大きい gcd を含む例で、バイト列と復元係数を独立に照合した。非空 q=0の拒否も確認した。元 worker から AST で取り出した `recover_committed` だけを試し、未完 manifest と hash不一致 tail が確定 prefix から除外されること、二世代の確定バイト列が一致することを確認した。選択 VALUE 関数や元の係数式の演算は呼び出していない。

数値所有量は、グラフ・numeric key・needed/order・thin segment index・runtime/controllerを再計数した。旧4096 control枠を除く baseline107,657に、同時に存在し得る二つの2000数値 receipt、129個の full5-key set、gate/controlと一時領域を加えた保守上界は **113,478< META118,000** である。各 VALUE map の項は係数と全5座標を含め7数値として数える。子 map は nodeごとに解放し、全1884の係数表を同時保持しない。

phase の範囲は companion により明確になった。通常成功時の最大は node141 の input decodeで138,535。child decodeの追加2048を含む最大は126,982、通常codecは125,192である。input decodeを保持した停止処理が partial codecを作る際には138,658となる。FULL F11 の anchor check時には120項の子 mapもまだ保持されているため、anchor停止時は128,892と数える。すべて元の hard150,000以内であり、gate/workerを上書きせず、この metadata companion を追加した。停止時も、正しくハッシュ検査された manifest/stream の確定 prefix と canonical temp cellを保持する。

整数高さは source M≤824ビット、q_iの分子≤28・分母≤164ビットを使う。G≤Mなので、codec value hashの素朴な一時積1676ビット、d_iG復元1648ビット、stored content/q_iの交差積1016ビットはいずれも1921<2048に収まる。入力係数と各演算の全 partial bound を source段階で確認し、実コードは M_i と1024bit commit guardを照合する。

一世代の cell最大9,268,352バイト、同一の canonical temp files と immutable streamの二世代を両方数える。manifest、root/control、停止時 partial/tail reserveを加えた generated上界は23,968,000バイト、固定入力・receipt reserveを含むfresh上界は26,155,445<32MiB。二世代の片方を消したり、圧縮後だけを数えたりしない。50秒の演算 guard、55秒 signal、60秒 hard、2GBを維持する。

独立 SOURCE 監査自体の owned上界は92,290、実行0.157秒だった。新しい係数応答、parameter point、prime、matrix/RHS更新はない。source支持129個という結果を実非零と解釈しない。実非零・階数・pivotの採用には、今後の ROOT許可実行と別の actual peer が必要である。既存worker/gate/method/reportのバイト列は変更していない。
