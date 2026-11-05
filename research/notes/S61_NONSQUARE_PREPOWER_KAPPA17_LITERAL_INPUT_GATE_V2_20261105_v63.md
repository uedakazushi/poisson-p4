# κ17 source 入力 guard の V2

2026-11-05。初回の source-only assertion は未使用の pair JSON を期待して停止した。
polynomial 値の評価前であり、失敗 code と receipt を別の旧ディレクトリへ保存した。
V2 は同じ20 atomの原num/denをコピーする。参照元はQ2とcheckpointの二つだけで、
係数record数40+101604、切出し2858、全metadata reserve35000から、
source準備全phaseは139502 slot。後続workerはこの小入力だけを読み、
元checkpointを再loadしない。convolution/publicationは75410 slot、整数temporary548bit。
旧SOURCE_GATEの75390はconvolution phaseのみで、旧111610はmatrix25x32のみを
数えた準備値だった。元資料は変更しない。

κ17や四source値は未評価。原20 atomを文字通りコピーするだけであり、scalar workerの
別GOが必要。零/非零どちらも最終一対八のrankや全FULLU57を証明しない。
