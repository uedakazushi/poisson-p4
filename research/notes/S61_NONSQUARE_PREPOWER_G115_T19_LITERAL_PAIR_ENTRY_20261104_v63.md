# 固定G115/target19 bordered六minorの一回計算

2026-11-04。先に凍結したsource-only見積もり108513係数枠・15.32MBに基づく一回の
QQ係数計算である。他行・他列、D5の展開、点・素数の追加は行わない。

最初に既存payloadからG20のanchor22が156項非零、g=N[6,1]が26項非零と確認する。
pair-pure六列のrank一と既存Saux非零により、D5_G20も一般点で非零となる。
これらは選んだ追加chartの非空性を示すもので、FULLUの大域単元性を示さない。

P4四行、literal G20、Ksmall=g G115/c10−h G85/c6を、元target順
[22,11,13,25,15,19]で並べた六minor D6を扱う。行の既知a因子は
[0,1,1,0,0,0]なので、正規化workerの出力にa二乗を掛けると元D6となる。
B poleのclearは全行で零、QQ分母とsourceで既知のa-Laurent offsetは別々に保持する。

1325の必要source演算を最後の使用で解放し、P4の十五の四minorとbottom二minorを
一つの整数格子へ融合加算する。十五個の共通QQ分母scaleが整数であること、
実contentを用いたscaleも整数であること、全source支持上界を実行時に検査する。
小さな独立QQ多項式積で整数格子実装を検査してから実sourceを読む。
native五十秒guardに達した場合は完了したLaplace項数と実部分和のみを保存し、
全minorの値とは主張しない。同じ計算の自動再試行をしない。

完了して非零でも、結論は追加D(g D5_G20)上の指定一成分とそのsource断片である。
全四残余行rank、全FULLU三次57、quartic239や尽くしは未決のままである。
full57 liftの独立peerを入力へ固定し、元八tailのうち既消費一列以外の七列を保持する。
