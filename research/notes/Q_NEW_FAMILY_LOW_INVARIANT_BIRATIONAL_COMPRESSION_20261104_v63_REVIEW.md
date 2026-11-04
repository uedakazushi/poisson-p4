# 五座標への再圧縮とB回収：逆代入による独立監査

2026-11-04。対象は frozen
`Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63`。
既受理のH約分後の七写像を使う。producerやその順方向変換を呼ばず、別の
Fraction疎算術により **新しい七多項式を旧座標へ逆代入**して全係数を比較する。
新しい点、source26、方程式探索、GBは行わない。

既知involutionの固定体は既受理の
`QQ(a=u/X,Y=X²,k=c/X,W,B)` である。

    xi=40k/(Ba),  r=a²Y/1600,  s=48aY,  q=3072W−1
    a=76800r/s,  Y=s²/(3686400r),  k=1920xiBr/s,
    W=(q+1)/3072

について、両向きの全合成を厳密に確認する。前後三変数とBの合成はLaurent
monomialの整数指数・有理係数で、q/Wは一次式として比較する。従って
`QQ(xi,r,s,q,B)` は **既知involution固定source体の全体**である。
actual NEW target の関数体と等しいこと、有限指数1は主張しない。

新しい多項式を `Ni(r,s,q,B)` とし、写像は
`xi^ci Ni / h^hi`（整数分母はNiのQQ係数へ含める）と書く。
逆代入 `r=a²Y/1600,s=48aY,q=3072W−1` を行う。
`h=(25r)H=a²Y H/64` と `xi=40k/(Ba)` を使い、非零分母を払った
**多項式等式**として、保存済み旧商qiと比較する。全七係数が一致する。
作者の W=(q+1)/3072 による順方向二項展開とは逆向きの検査である。

項数は `5,6,29,33,31,30,43`、残るh分母は `0,1,1,1,1,1,2`。
各整数分母、xi指数、全multidegreeも確認する。

    h=s²−25rq²,
    F=100r(r−1)q²−4rs+s²+4r,
    F4=(2/405)xi²F

および p0 の六項分子の因子表示を独立に乗算確認する。さらに
`p0+4p4=F4²/400` を同じ七写像のQQ多項式として検査する。
従って Q0=p0/F4² に対し

    Q0 = D B/(hF²),
    D=2359296000r⁴(r−1)q²(−s+5q+1),
    B=Q0 hF²/D.

比例定数も厳密な有理数として比較する。Dはr,s,qの非零多項式なので、これは
非空なdense rational chartでのB回収である。必要な逆元は表示しており、
`r*s*B*h*F*D≠0` と元source domain（scalar非零を含む）に制限する。
既知固定source体ではBをQ0へ置き換えられるが、r,s,q,xiをactual targetから
回復したわけではない。全非平方問題にこのlocalizationを適用することもしない。

この結果は既存七写像の小さい完全な表示と座標同型である。
NEW hypersurfaceの定義式・次数上界、全不変式環、target quotientの有理性や
birationalityを結論しない。専用reportは両向き座標合成、全七逆代入、因子・定数・
B回収の成否とhashを固定する。
