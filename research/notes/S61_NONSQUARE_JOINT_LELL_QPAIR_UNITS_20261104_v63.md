# 残る q0/q5 の58x136：既知単元だけによる一回の消去

入力は全Z/W補正と五つのoverlap単元を正確に処理した凍結済み58x136。
全183行列と全33genuine quarticは展開しない。
FULLU の a,d,F,B,C,J の monomial と非零QQ定数の積だけを pivot とする。
最初に a,d monomial の短い候補を調べ、fill見積り
`column nonzeros * row numerator terms` が小さいものを選ぶ。
候補がなければ残る既知FBCJ因子も正確な除算で確認する。
未知因子は反転しない。

全pivot値・既知単元factorization・全elimination係数を記録する。
選ばれた元58行と136列の minor は、選択順に並べれば pivot の積であり
既知単元である。各step完了時に全行列と全履歴をatomic checkpointへ保存する。
時間・項数・32MiBの制約で途中stepが止まれば、その未完成stepは数えず
直前の完成checkpointだけを有効状態とする。

一回60秒・2GB、intermediate項数上界150000、output各係数100000項。
追加の全域image pivotはkernel方向の追加ではないため、既知38・残り42は
変えない。full239の結論は与えず、残る全式を損失のないsourceに保持する。
