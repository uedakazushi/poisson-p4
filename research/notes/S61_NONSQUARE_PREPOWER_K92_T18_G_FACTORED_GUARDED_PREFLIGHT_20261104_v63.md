# gを外へ残すdet11と安全な疎prefix計算のsource gate

2026-11-04。前の二つの出力FAILは変更しない。gは旧N6 chart上で既に単元なので、
scalarの対象をdet11とし、実成分のc8 g/(Delta a二乗 N7)は因子のまま残す。
新しい多項式因数分解・gでの係数除法・parameter点の検査は行わない。

十H3係数のうち六個は既存のDelta/lambda値である。残る四個は元三対三H3行列から
支持・高さだけを見積もる。g H3minor=±h lambdaまたは±h120 lambdaという既存恒等式から、
支持函数の極値は積で加法的、QQ contentもGaussにより積で加法的（乗法的）である。
そのためgを実展開せずNewton区間と正QQ contentを得る。高さは除法で増えないと仮定せず、
元三minorの直接l1上界をこのcontentで正確に割った上界を使う。

完全出力の従来格子がまだ大きい場合に限り、次の安全な疎attemptを定義する。
固定した十H3項を全て足し、pure八minorは二層のsubset DPで計算する。積の多項式bufferを
作らず、integer addmulを一係数ずつ行う。固定QQ因子を各nodeへ割り当て、途中和でも全scaleが整数となる。
commit前に新cell数・係数bit・QQ因子bitを確認し、消えた零係数と不要DP stateを解放する。

全入力atom・DP state・重み・global和を同じ135000 coefficient cell poolに数え、
15000作業枠との合計を150000とする。active polynomialは1024、各整数係数1024bit、
共通QQ因子4096bitを上限とし、越える直前で止める。これは完了を保証する仮定ではない。
native clockは45秒で止め、残りを正確なprefix保存に使う。

停止時は全active polynomial、common QQ因子、row/maskと畳み込みcursorを一度だけ
binary保存する。係数一件の上界144byteと共通因子/headerを数え、入力・metadata・log・
manifestも含めて32MiB以下とする。古い巨大checkpointと新temporaryを同時保持せず、
係数をJSON/base64へ重複保存しない。未完prefixの非零をdet11の非零とは言わない。

このgateはsource/資源設計だけであり、scalar実行はしていない。
全十項・a lambda419・tail2効果・全補集合を保持し、全FULLU57へ推論しない。
