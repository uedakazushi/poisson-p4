# exact37残余のindexed row-stream representation：source-only gate

2026-11-04。入力は凍結された199×46、147370項の完全な論理行列である。
新しい算術pivot、係数の積・差、gcd、級数、点・素数は実行しない。
旧all-resident guard missも全過去の結果も変更しない。

## 新しいguardの意味

150000は引き続き同時に所有する係数とscratchの上限である。論理行列の全展開項数は
別に報告する。RAMには一つのcurrent row、一つの固定pivot row、一つのentry更新の
作業領域と小さい候補metadataしか置かず、他のliteral係数はimmutable fileに全て残す。
これは係数を捨てる操作や上限の拡大ではなく、保存representationの変更である。

## 全literal入力のindex

元JSONを64KiB raw chunkで読み、同時には一つのcoefficient-cellだけをdecodeする。
全5358cellのbyte offset、length、payload hash、元行・列番号と項数を記録し、
総項数147370と元file SHAを確認する。全zero rowとactive列IDも保持する。
入力全体をQP配列としてmaterializeしない。

row readerは必要なpayloadだけをseek/readし、hashを確認して一つの行を復元する。
全行のterm数・byte数はmetadataで持ち、係数値は持たない。単項式候補のQQ係数
metadataも保守的にowned予算へ算入する。file index自体もRAW保存量へ課金する。

## first候補の完全なsource/height予測

選択はmetadataだけでpivot row項数、target read volume、元行・列順を最小にする。
同じliteral parameter単項式を使い、非単項式field inverseは許さない。
全target行を一つずつ読み、全entryについてsupport和集合、numeric分母LCM、
primitive-content桁数不等式、実際のLaurent指数から予測する。

current rowの予測所有項数をLrow、固定pivot rowをRとし、候補metadata分も加える。
既検証のexact Laurent product/sub/back-add kernelについて

    Lrow + R + metadata + 2k + 3b + 96

を各更新前に確認する。rowを解放すれば次rowへこの小さい所有状態から移れる。
first gateでは新しいfused arithmeticを仮定せず、既検証kernelだけで予算を証明する。
全論理matrixのsupport和集合上界と、旧all-resident予測も別欄に保存する。

## 全係数を保持するgeneration

後に算術workerが承認された場合、全row更新とback-addを検証して変更rowを
新しいimmutable fileへ保存する。非変更rowはhash-boundな元index参照を再利用する。
正規化pivot rowとsource factor参照、元行・列順をreceiptに残す。
全target rowが完成してからだけ、新しい完全manifestをcommitする。
途中停止ならold manifestが有効なままで、完成済みpending row・identityを保存し、
新しいpivot完了とは主張しない。

一時generationについて本gateはdisk削除を必要としない。parsed row/cellはRAMから
解放するが、jobが作った全row versionをRAW ledgerに残して課金する。
将来disk versionを消すには、残したinputとexact operation pathから全係数を再構築
できる別の検証が必要であり、externalと呼び替えて32MiBを回避してはいけない。

RAWにはcurrent frozen inputを一回、literal index、新しい全変更row、正規化row、
全receipt/preflight、manifest/cursor reserveを含める。負指数を捨てず、全五parameter、
全alpha/beta制約と全論理係数を保つ。密inverseも作らない。

## IOと時間の限定

source gateはhard60秒・2GBで測定し、読み込むpayload bytesと将来のread/write上界、
係数pair数も保存する。byte数はOSのthroughput保証ではない。
将来の算術workerにはdeadlineとsoft stop、上のatomic manifestが別に必要である。
このsource-only gateだけでは算術GOを得たとはしない。
