# 全Qa上のCA/CB/N proper-source微分：QQ[a] Smith証明書

2026年11月3日。これは**source differentialの像**の階数であり、target Zariski tangentの次元ではない。

物理座標(h,x,y,z,w)=(old0,old1,old4,old2,old3)で

    H=(4/3)hxw−hz²/3−x²z,
    Qa=raw_CA(a*h²,−2/9,H)=raw_CB(a*h²,−2h²/9,H)

とする。全150 canonical係数の像微分を次の全fiber方向と全GL25方向で作った。

* CA：binary fの3方向、g h²の1方向、Hの16方向。
* CB：pure f=h²、pure g=h²、paired(f,g)=(a*hx,−2hx/9)、(a*x²,−2x²/9)、H16。
* N：f=h²,hx、g h²、H16。

weight(h,x,y,z,w)=(0,1,4,2,3)で全行列をblock分解し、各blockのQQ[a] Smith形と
unimodularな左右変換を保存した。全積恒等式と左右determinantの非零定数性を検査した。
従って以下は全複素aについて成り立つ。

|source|a≠0のaffine image rank|a=0のrank|非定数Smith因子|
|---|---:|---:|---|
|CA|28|28|なし|
|CB|28|27|a/9が一つ|
|N|27|27|なし|

a=1とa=2/9の間で、これらのproper-source像微分は変化しない。
actual target成分のcone次元29,29,28よりそれぞれ一つ低いので、source微分だけからsmoothnessを結論してはいない。
a=0は物理的rank2点であり、S61のsymplectic-value正規化とは別に保持する。

* `code/q_CA_CB_N_source_smith_20261103_v63.py`
* `reports/Q_CA_CB_N_SOURCE_SMITH_20261103_v63/{report.json,certificate.json.gz}`

全QQ[a]計算1.84秒。全150行列・weight分解・全左右Smith変換を保存した。
