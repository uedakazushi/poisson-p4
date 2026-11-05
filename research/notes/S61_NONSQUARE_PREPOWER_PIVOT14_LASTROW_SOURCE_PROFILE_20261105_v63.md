# 一対八の残余と未処理 tail の実 source profile

2026-11-05。係数の新評価は行わず、accepted Q2/L completion の有限 DAG から
current15 全行と残る八列 [1,3,4,6,7,8,9,17] の依存閉包を記録する。
PIVOT14 の非空性は別の独立非零確認を待つ条件付き source である。

実 source の rawpair 四行は current15 [2,4,5,8]。
既証の rank 一は固定八列 [22,11,13,25,15,19,14,18] に限られる。
今回残る八列はいずれもその外であり、三つの pair 行を無補正で捨てて
M11 と G134 だけの十二角 border にすることは正当化されない。
Q9 の全 tail、消費済み tail2、および三つの pair-row correction を保持する。

current15 全十五行の closure は target17 で 159 node、64 atom、原 numerator
10406 項であり、最大 atom は 394 項。rawpair 四行だけなら49 node/20 atom/
2838 項、最後の G134 一行だけなら19 node/8 atom/1254 項となる。
他の七列は全十五行の原 numerator が約3.36〜3.52万項である。
これらは原 atom の件数であり、DAG 出力の展開項数・行列式の仕事量上界ではない。

全 denominator を QQ 非零定数と a の単項式、または B の定数倍として
文字通り確認する。未知の chart 因子を atom に導入しない。
各 closure、root、atom reference、項数・指数 box・bit 高さを保存する。
後続の実 entry は PIVOT14 の cleared root s*d-c*b を使う有限作用素である。
少ない原 atom だけを理由に展開 worker を承認したとはしない。

次の実計算には、この cleared root の全既存逆元を tracked parameter 因子で払う
具体的な和・積 schedule、全 lifetime/整数高さ/RAW 上界を別に示す必要がある。
新しい係数、minor、rank、点、素数は本 packet では一切計算しない。
全 FULLU57/239 と chart の補集合は開いたままである。
