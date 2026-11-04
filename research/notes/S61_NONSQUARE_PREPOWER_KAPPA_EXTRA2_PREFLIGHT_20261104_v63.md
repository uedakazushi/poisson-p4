# 射影source κ・追加tail列2の厳密な高さと生存支持

2026-11-04。SOURCE-ONLY。一つの列2だけを扱い、parameter係数は新たに計算しない。
四sourceはcurrent15 `[2,4,5,8]`。H二列はFULLU既知単位による正規化済み、anchorは
既存column22。符号付き四cofactor κと全四行の合成をそのまま保持する。

列2のsource原子にはliteral QQ*B分母がある。全分母を逐語検査し、各中間RFを
`numerator/B^b`として型付けする。加減算では共通B冪へ上げるために必要なB多項式を
明示nodeとして掛け、最後に全四行の共通B冪を掛けた列をκと合成する。未知分母や
parameter点での非零性を仮定しない。native候補へ渡す全nodeはa-Laurent/QQだけになる。

各原子の全monomial、QQ分母、共通分母を払った整数l1上界を読む。積では次数・指数boxを
加え、和ではそれらの和集合を取り、全degree区間内の整数格子点を数える。l1=0又は同一
node差だけを零とし、未知cancellationは仮定しない。各係数の整数height、QQ分母、
source/cofactor/outputを保存した場合のserialization上界を算出する。

全nodeのlast-use、preloaded literal atoms、output-sized scratch、15000 work reserveを
含む150000項の実装見積を保存する。四cofactorは一つずつ計算して使い、全四行の合成を
減らさない。全8追加tailの保持を元sourceへ明記し、今回評価対象はそのうちの一列のみ。

対象は実際に保存された**projected prePOWER15×22**である。ここでのpureは表示された
heavy5の消去を意味する。元raw H7のpure sourceを主張するには、既存current15→old25と
eta単位によるHZ/HW/Cdelta補正を別に合成する必要があり、その未展開gateを保持する。
係数が零でも非零でもfull cubic57や有限bad locusを結論しない。
