# CA V1でのgeneric V cusp-order：実際のodd制約を除く有限計算

2026-11-03。これは実行中の計算checkpointであり、cusp次数や全局所型を証明済みとはしない。
CBのQQ(V)計算はまだ開始していない。

## raw target sliceを使わない

目的はgeneric V境界でのactual odd idealから八つのregular generatorを除いた後、
残る一変数nuの最初の非零次数をp³,p⁵,p⁷で調べることである。
raw28 target座標Y=0をそのまま固定すると、odd regular constraintの高次補正を落とすため不適切である。

sourceの29変数をG20＋fixedfiber補9で取る。
target微分のkernelをkappa、他28方向をrとし、source(p,r)=kappa*p+S*rと書く。
20個のactual orbit座標を固定し、29個のactual quadratic-base出力と選択八odd出力が
p→−pでevenになるようpaired sourceを構成する。

source kernel coordinateのodd部分をepsilonだけに固定する。
偶数次数2dのsource係数29個と、次の奇数次数2d+1の非kernel係数28個を未知数とする。
課す条件は

    20 orbit出力のeven係数 +29 base出力のodd係数 +8 odd出力のodd係数

の計57本で、毎回同じ57×57有理行列を使う。
このblockの可逆性と別の28次minorにより、選んだ八元がactual odd regular constraintsであることを確認する。
残り全150出力のodd係数を次数3,5,7で検査する。
pairのbase出力がevenでsource pのodd部分がepsilonなので、tracefree t=epsilon、nu=epsilon²となる。
したがって残りのodd次数2k+1が、八regular odd元で割ったnu-idealの次数kを測る。
全R-baseのfixedと八odd制約は区別する。

## 最初のbounded runは失敗した

`code/ambient_CA_V1_odd_order_20261103_v63.py` は60秒上限で停止した。
maxRSSは155900KBで、cusp-orderの結果はない。
**setup57は保存されていなかった。**
凍結failure reportの一般的なnext_step文にある「saved setup57を再利用」は誤りであり、実在しないものを使ってはならない。
元source・失敗reportは改変せず、本追記で訂正する。

## 保存できたcheap numeric frame

`code/ambient_CA_V1_odd_order_checkpoint_20261103_v63.py` を1.39秒で実行し、
次を`reports/AMBIENT_CA_V1_ODD_ORDER_CHECKPOINT_20261103_v63/checkpoint.json.gz`へ保存した。

* V1の全physical150係数とfixedCA20基底。
* 選択G20のGL generator番号と補fiber9の基底番号。
* 全normalizer source29 frame、actual target微分150×29、kernel29、complement28。

orbit rank20、normalizer source29 frameは可逆、target source微分rank28である。
これら以降の巨大SymPy Hessian/57solveを再実行していない。

## 最適化した独立作業

親agentの指示により、v62 reviewerが疎な150係数Lie作用から数値mixed Hessian H_kappa,jを担当する。
新runner `code/ambient_CA_V1_odd_order_v2_20261103_v63.py` はこのhelperと保存frameを読み、
source Gをordered formal exp(s_i E_i)で評価する。
旧失敗実装のordered(I+s_i E_i)とは対角方向の二次以降が違うため、
helperもrunnerも必ず新しいformal-exp規約へ合わせる。Jとkernelは変わらない。

runnerは二次source係数とH_kappa,kappa/2を全150で比較し、
57単元blockの補正も実source再評価で確認する。
次数3,5,7の検査はunder60s/2.2GBの新しいbounded runであり、失敗したらhonest failureとして保存する。
成功してもgeneric Vの全局所型には、dense orbitとformal group sliceの別の適用証明が必要である。
