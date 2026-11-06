# 元134式と固定etaのnative支持判定：入力と到達範囲

2026-11-06 / v64-classification-first。Gの計算入力を固定した。以下は保存された一つの有理一次jetを法32003へ移した診断であり、全K16上のeta所属や一般五点尽くしの証明ではない。

## 数学的な入力

全134式の元61変数エンコーディングへ、保存T01点の一次jetを代入する。45個の二次座標だけは未知数のまま残し、q=q_T01+deltaqへ正確に中心化する。既存のpulledfiber helperをproducerとして再利用せず、各元単項式をQQで展開して構成した。別実装の全134式との一致をsource peerが確認した。

採用するnormalはglobal J=[0,1,2]のframeであり、

    eta=n0+n34
       =deltaq3+deltaq37-(Wrow3+Wrow37)*deltaqfirst3.

疎K5 typed模型の同じ番号のnormalへ置き換えない。135本目として1-u*etaを加える。全134式を保持し、元の二次項を捨てない。基礎体は今回GF(32003)であり、全QQ係数の分母がこの素数で単元であることを確認した。

## G6によるexact前処理

既存の`S61_UNIVERSAL_WEIGHT6_20261031.json`にある125項のcofactorをc_jiとすると

    G6 = sum_(j,i) c_ji * u_qj * F_i,   i=0,...,69

は元ideal内の恒等式である。新しいPoisson条件を追加するわけではない。

今回、元134式から全9257個の普遍QQ係数を作り直し、保存したG6と全項一致させた。次に中心化し、790項の二次以下の式を得た。

    G6(deltaq)=sum_i K_i(deltaq)*F_i(deltaq),
    K_i=sum_j c_ji*(deltaq_j+q_T01,j)

もQQで全係数を確認した。source peerは別に普遍9257項・固定790項・K_iを再構築し、元135本の不変性とSingular入力の全136本のGF係数を直接parseして一致させた。

`classification_native_support.py`はprepare/runを持つ一つの再利用runnerである。--add-g6は上のexact ideal consequenceを先頭へ置き、--algorithmでstd/slimgbを選ぶ。いずれも全134式と同じetaを保持する。

## 実行結果の限界

* 元135本のstdは300 CPU秒・2GiB・315 wall秒のwindow内に終了しなかった。最後までbasis完成markerは無く、unit・nonunitのどちらも得ていない。終了codeが0でも、Singularの終了処理によるものであり成功とはしない。
* G6を加えた136本のslimgbは、600 CPU秒・4GiB・615 wall秒のwindowを設けたが、実際には約33.0 CPU秒・33.0 wall秒で`no more memory; halt14`となった。peak RSSは4,172,232KiB。basis・unit・nonunitはいずれも未決である。

* 同じ136式・同じSingular scriptのbyteを保った最終12GiB評価では、600 CPU秒・615 wall秒を設けた。起動直前にcgroup上限16GiB・使用817,397,760 byteを読み、12GiBを差し引いても3,477,569,536 byteの余裕があり、要求した2GiB reserveを満たすことを確認した。しかし約64.8 CPU秒・64.8 wall秒で再び`no more memory; halt14`となった。peak RSSは12,560,924KiB。basis・unit・nonunitは未決であり、追加の自動再試行はしない。

これらの失敗は数学的な非所属や残余存在の証明ではない。slimgbについては時間だけを延ばせば解けるという結果ではなく、現環境で採用した方法の実際のメモリ制約を示した。異なるbeta・素数・chartを試してはいない。

小さい結果reportは`CLASSIFICATION_NATIVE_STD_FIXED_GF_20261106.json`、`CLASSIFICATION_NATIVE_G6_SLIMGB_FIXED_GF_20261106.json`、`CLASSIFICATION_NATIVE_G6_SLIMGB12_FIXED_GF_20261106.json`。最終資源契約は`research/classification/native_slimgb12_job.json`。大量のplain入力・log・厳密source hashはignored workに保持し、runnerから再生成できる。native環境の版と導入記録は`CLASSIFICATION_NATIVE_ENV_20261106.json`。


## frozen betaの既知例外像について

別担当のexact QQ評価と逆順整数和による独立確認は

    G28(beta_frozen)=671334843179656461002307382635849/9039207968 != 0

を与えた。同じ法32003では分母が単元で値は25057≠0である。従って今回の固定betaはQQでもGFでも既知NEW25の閉像から外れている。この確認は他の48既知proper像や全未知例外像を排除しない。

証拠は`CLASSIFICATION_FROZEN_T01_BETA_G28_20261106.json`と`CLASSIFICATION_FROZEN_T01_G28_FINAL_SUM_REVIEW_20261106.json`。このpremiseを追加してから最終12GiB評価を許可した。今回のnative jobでこのG28値を再計算したわけではない。

## 証明書の扱い

basisが1になった場合だけ、同じnative CPU window内でliftを試みる。完全なcofactor exportが得られたときは、別のPython疎積算法で1への掛け戻しを確認してからplain JSON証明書を保存する。G6のcofactorは可能なら元135本へ展開してもう一度確認する。確認予算を超える場合でも、136本の掛け戻しと上の明示G6恒等式の組はexact証明書となる。

今回はそのようなunit証明書を得ておらず、数学的証明書の空ファイルも作成していない。仮に固定GF fiberのunitを得ても、全K16上のunitや標数零の一般空性へは自動的に移らない。

次の構造的な候補は、principal eta上で全cross solveから三つのzを先に消し、全compatibilityとcleared Fを残したnormal-only表示を作ることである。これは今回まだ実装していない。etaを単に係数体の超越パラメータとして扱うと、etaが非零定数である成分を失い得るため、その変更で代用できない。

## 再現入口

新しい出力先を用いる。実行にはlocal Singularとそのoptimized module環境が必要である。

```sh
.venv/bin/python code/classification_native_support.py prepare \
  --output work/new-native-source --field gf --algorithm slimgb --add-g6
```

VALUEはsource、全hash、有限job契約を確認してから、同runnerのrunへ新しい出力先・source hash・resource limitを渡す。過去の失敗を理由に、同一jobを自動反復する仕組みは設けていない。
