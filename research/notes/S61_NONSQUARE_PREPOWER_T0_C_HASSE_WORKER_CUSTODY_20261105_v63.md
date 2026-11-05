# 旧T0だけのC-Hasse workerと入力custody：未実行

2026-11-05。既存KAPPA_EXTRA2 root239のT0を使い、C-valuationの固定workerを準備した。
元の7208係数は全QQ source peerに固定済みであり、新しいtail・点を選んでいない。
compileとhelpのみ実行し、C冪・Hasse層・valuationは未実行である。
fresh ROOT/source-code全体監査と一回GOを実行前に要する。

root239.polyのheader4はa-offset。bodyは既にU0=T0/a⁴のa指数0..20であり、
**将来workerのshiftは0**とする。変数aそのものをa−4へ置換する操作でもない。
bodyの7208 monomial順を保持し、既存共通D197でscalarだけをI=n*(D/d)へ直す。
各d|Dとscale*d=Dを検査するので、大きなQQ cross-productを要さない。
integer最大320 bits、l1は326 bits。

C=100a+18e²−9e−2、a0=(2+9e−18e²)/100。
k0..20の最初の非零だけを求める。一個のaccumulatorでforward／逆順減算零／forward再現を行う。
input7208、root cache441、e≤80,d≤18の1539 slotsと全metadata16000により、
coefとpacked key双方を数えて34376 cellsを維持する。
input_l1*2²⁰*100²⁰からHHの全partialは481 bits以内、hard518は上げない。

全conversion/scaled/root/products上界4117089。
45秒soft/60秒hard・2GB・RAW8MB。inputvectorは7208をreserveし、L1もcommit前に検査する。
prefix T0CHP001は全input/cache/accumulatorと全cursorを一ファイルに保存する。
partial prefixからvaluationを主張せず、自動再開しない。
出力I/O失敗に第二のfull prefixを重ねない。

この準備は他因子、GCD、他tailやbad locusへの結論を含まない。
既存FULLU・D12の十分条件と全未処理補集合をそのまま保持する。
