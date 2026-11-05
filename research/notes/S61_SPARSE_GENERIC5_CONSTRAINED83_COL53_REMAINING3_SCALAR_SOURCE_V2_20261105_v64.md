# LIVE07の既存支持検査と残3成分のscalar/BOX上界

これはSOURCE支持計算の再開でもVALUEでもない。既存の第7境界snapshot22map・44555個の5byte keyを独自readerで一mapずつ読み、重複・QQ/M・個数・全BOX・Hcountsとaccumulator canonical hashを検査する。元snapshotや旧guard155866を変更しない。

そのmetadataと既存leafだけから、残るraw成分8,9,10（A11,A44,A52）の93項についてQQ scalar、正M、H区間だけを伝播する。新しいK5多項式の係数もMinkowski支持も評価しない。materializeする一つの支持は、同じ150k SOURCE guard・base最低25kを通る場合に限り125000以下として上界化する。未materialize dotのcardはこの上限で切らず、積のpotential pair数を使う。Mは各出力への寄与個数min(nA,nB)と正の入力Mで抑える。

この上界は、immutable support＋offset viewによるSOURCEコピーと、唯一所有accumulatorへのSOURCE Minkowski union融合を検討するためのもの。QQの共通scalarに対するwa/wbはLCM/gcdだけから計算し、qNew*wa=qOld、qNew*wb=qA*qBを保つ。SOURCE dotを未materializeにする場合、そのcard/Hcountsはupperと明示し、最終accumulatorのcard/Hcounts/hashはexactと区別する必要がある。

新SOURCE workerの所有数・全RAW・storage往復・code-review gateは別であり、このreportから計算再開やVALUEを許可しない。全167×14、87追加RHS、66制約RHS、P290一回、F9/F10と全補集合を保持する。

V1独自snapshot readerは空supportのBOXを[191,−64]の初期値のまま比較して停止した。SOURCE protocolの空BOXは[0,0]である。V1code/failureは保存し、V2で空集合だけを正規化する。SOURCE profileやVALUEは実行しておらず、snapshot入力・式・予算を変更していない。
