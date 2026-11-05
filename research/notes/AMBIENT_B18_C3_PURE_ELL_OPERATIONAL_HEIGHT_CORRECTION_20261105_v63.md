# Pure-ell cache replay の演算途中 height guard

元 SOURCE の2048-bit 継承は、Q20で保存された scalar DAG coefficient に関するものだった。小さい再出力 worker でその係数を得る途中の Fraction 加算・畳み込み・単元逆の診断積まで同じ bound と解釈するには、別の operational check が必要である。この点は値 job の前に独立監査で見つかり、値は一切計算していない。

frozen V1 と custody はそのまま保存し、別名 V2 を作る。各 Fraction x=n/d,y=m/e に対し、積の numerator/denominator の bit 上界を max(bits(n)+bits(m),bits(d)+bits(e))、和を max(max(bits(n)+bits(e),bits(m)+bits(d))+1,bits(d)+bits(e)) として、演算前に2048以内か判定する。零の項はそのまま省き、負号は bit を増やさない。

全3係数の convolution、六つの原点1単元逆の式、その inverse product check、physical delta の引き算もこの checked helper を通す。入力の分子・分母は Fraction を作る前に確認する。上界を超えた場合は未完の requested root だけを checkpoint に残し、当該演算の前で停止する。cap を上げたり値を0にしたりしない。

source、入力、式、対象164 root、全105式、gauge/IFTを解かない範囲、全容量限界は変わらない。controller の増分は既存10000 reserve 内である。pure helper fixture と import-only --help のみを確認し、実 source readback はまだ実行していない。
