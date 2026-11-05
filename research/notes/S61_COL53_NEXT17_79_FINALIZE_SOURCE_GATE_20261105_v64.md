# 保存済み SOURCE prefix の read-only 確定ゲート（未実行）

旧157050停止は変更しない。新workerは、凍結した全2998ノードprefixを一回読み、補助 box5 と n_upper のみを resident view から除く。元prefixと全旧履歴を保存するため、正上界の根拠は失われない。旧コンパイラのQQ/GCD/係数側演算は一つも再実行しない。

stream census の86,903から、非null BOX 2320個×10および n_upper 2998個を除き、resident graph の上界は60,705となる。cold parseは115,999、依存map・stack・Root文書・薄いpayload alias等を含む確定処理は125,689、固定上限は126,000である。全生成RAW4MiB、fresh全読込・出力RAW6,938,295バイトで、150k／32MiB／55秒signal・60秒hard／2GBを維持する。

依存索引は共通ノードの最大相対s次数を取り、正の先頭H,Lとs範囲による枝だけを落とす。実際の頂点零・面零で先頭値を動かさない。serialized active mapの値listは元mapとの同一aliasであり、係数・指数配列を複製しない。mode map・serialized key/indexに8N、stack等に12k、追加custody等に4kを見込む。

旧raw53の入力は、Source前処理で保存済みの同一payload SHAを照合して薄いpackへ参照をつなぎ替える。元14MB行列の履歴と元参照を保持し、新jobではそのファイルを読み直さない。これはbyte同値aliasであって新しい係数値ではない。

workerの --finalize はゲートSHAとRootレビューSHAを必須とする。現時点はゲート準備のみである。確定後も新しい射影支持集合・符号付き係数は未計算で、別ゲートを要する。全166×13・87追加RHS・66制約・P290一回・F9/F10/F11を含む全補集合は保持する。
