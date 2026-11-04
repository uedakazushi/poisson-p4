# 残る四つの線形Q例外成分の全有限代数検査

factor index1..4だけを一つのbudgetで調べる。e−2/3、e+1/6、e+5/9、e+8/9。
各閉成分全体でQ(a,e0)を作り、a=0 clopen部分だけを既知単元として除く。
a有限代数とモニックHの基底a^i d^jを使い、全22元sourceとUを評価する。
これはgeneric点のsamplingではない。

全ideal columns Pi*basisをQQで作り、U^Nのmembershipと元22sourceのcofactorを
保存する。別の有限代数実装でsourceとcofactor dotを確認する。成功した成分では
FULLUにB-bad点が無い。失敗の場合は全行列を残し、globalなrank否定をしない。

以前のe+1成分と合わせ五線形成分が全て閉じても、e次数6,6,9,11の四成分と
旧有限h1/k1残余は未検査のまま。これらへの自動jobやGB・点・新素数は行わない。
