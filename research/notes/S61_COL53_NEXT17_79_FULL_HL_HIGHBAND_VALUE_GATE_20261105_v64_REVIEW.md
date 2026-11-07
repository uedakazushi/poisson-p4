# HIGH 帯だけを追加する VALUE の実行前独立検証

worker `6bba032b…`、gate `29c95901…`、custody `5ed3646c…` の SOURCE・コード・全体予算を独立に照合して PASS。今回、符号付き HIGH 係数は計算していない。既存の LOW 値1884個は、既受理の独立 local identity peer `2b0b2212…` を用い、元の stream/index を書き換えない。

各 node では SOURCE の固定 HL leading face を扱う。LOW はその node 自身の cap 以下、HIGH は cap より大きい部分であり、二つは互いに素で全 face を覆う。子の LOW と HIGH はそれぞれ独立に primitive 化されているので、保存 content を固定 qi で割った、そのバンド自身の正整数 G を掛け戻してから結合する。LOW と HIGH の G が同じだとは仮定しない。符号付き qi、異なる G、空バンドについて独立 codec fixture を確認した。空の場合は数学的 G=0、codec sentinel G=1 を区別し、qi による割算を行わない。

積では二つの子の完全な face を掛けてから親の HIGH 領域だけを取る。従って LOW×LOW から HIGH へ入る項も残る。和の枝は SOURCE の HL が一致し、固定 q が非零のものだけを使う。実際の leading coefficient の相殺や、以前の vertex 値が零だったことでは枝を切らない。VALUE は五つの指数と係数をすべて保持し、一項七数値の予算を使う。

独立に、元の literal leaf の枠と保存済み positive SOURCE 子 support から、全1884個の完全 positive support を再構成した。494,806組の Minkowski 積を含み、全 node で保存 support と一致する。LOW の positive support が同じ cap による切断であることも全 node で照合した。全 HIGH positive slot は191,169、root の完全 positive support は1012、HIGH は883。これらは係数の非零数ではなく SOURCE 上界である。

各 COPY・ADD・MULTIPLY の固定 qi 正規化比が整数であること、全 face の Cauchy 係数上界と途中和の絶対値が node の保存 M 以下であることを再計算した。最大 M は824 bits。整数 rescale weight の実際の最大は31 bits、node2973の −1134000000 である。以前のレビューの例示28 bitsをここへ流用せず、保守的な積上界 `2*824+31=1679 <1921 <2048` を使う。この違いは既存の1921 temporary guardを変更しない。最初の独立 SOURCE draft の上界 assertion とコードは履歴に残した。

metadata は未使用の active-HL 補助表だけをメモリから退役させる。元の graph ファイル、全 scalar node、SOURCE lead、実 LOW は不変。resident graph は62,623、mode keyとorderが各1884、combined index は key、LOW三欄、HIGH三欄、二つの hash を含む9Nである。runtime25,000、現 node の制御4096と合わせた census は112,443。二つの将来レビュー文書4000と余白1557を含めて META118,000に収まる。plan と旧 manifest は一行ずつ読み、全 plan または第二の旧 index を別配列として持たない。

全1884状態の独立再計算で、通常時の最大は143,570、子 decode 時に root2997で生じる。leaf decodeは138,542、子を解放した後の codecは130,926、anchorは128,629。途中停止時は leaf と codec が同時に残る138,668が最大で、anchor停止は129,142。子は通常の codec/anchor 前、および例外時の不完全 accumulator codec 前の両方で解放する。すべて hard150,000以下である。50秒の数学処理期限、55秒signal、60秒hard deadline、2GB ASを保つ。

新しい HIGH 保存は canonical cell と stream の二世代である。codec 一世代の上界10,676,690 bytes、二世代・1884 manifest行・停止制御を合わせて26,784,676 bytes。実 READ3,953,454 bytesと将来の二つの64KiB文書・制御余白を含む fresh上界31,647,203 bytesは32MiB以下。manifest一行の最悪値1749 bytesも固定2048以下である。元 LOW 二世代を追加の新規生成として数えず、旧 stream/indexは読み取り専用として全体SHAを再確認する。

独立 SOURCE 検証の所有上界は145,045数値、39.69秒で完了した。これは未選択 HIGH 値の評価ではなく、support・正の上界・入力枠・コード・無関係な codec fixture の確認である。実 VALUE は親の別の ONE GO を要する。将来 root face が零でも whole minor の零や rank を主張せず、非零なら保存済み full five-coordinate witness に限って結論する。新しい matrix/RHS 更新、パラメータ点、素数、元 source 係数の solve はない。
