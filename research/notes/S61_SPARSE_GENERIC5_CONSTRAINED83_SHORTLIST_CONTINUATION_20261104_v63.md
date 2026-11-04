# 30段後のsafe保存候補から、shortlistで進めた6単元対

2026-11-04。元から−t0を先に使った30段の厳密残余206×53・85459項を入力とする。
その30段も、別経路の旧24/29段も再実行・変更しない。同じreordered経路にだけ
今回の6段を足し、単元対は36個となる。

## source-only gateと変更した選択範囲

保存された未完了選択の125件には、完全にsafeと予測済みの候補が4個あった。
その中のwhole forecast最小値は元(12,67)=−32t0/405、90349項であった。
同じ凍結残余から全row updateを係数を作らずに再予測し、90349項を再現した。
符号付きLaurent shift、numeric分母LCM、primitive-content桁数不等式と実指数で
RAW上界も確認した。新しいRAW入力ファイルまで含めた上界は9248256bytes。
gateを保存して条件成立を確認した後にだけ、親が条件付き許可した一回を起動した。

以後は全literal定数／parameter単項式をcheap structural fill、項積majorant、
元行・列順で並べ、その最初の最大16候補だけをwhole forecastする。その中でsafeな
完成予測が最小の候補を選び、同点は元行・列順で決める。全候補での最小性は
主張しない。shortlistにsafeが無ければ停止し、未調査の候補をunsafeとは呼ばない。

## exact演算と予算

各全pivotの前にsupport和集合とL+2k+3b+64で所有係数とscratchを確認した。
全RAW入力も32MiBに算入し、receipt・exact残余・選択記録と次の未完了選択用reserve
まで含めて予測した。gzipは使わない。sourceを読んだPython配列はnative起動前に
解放し、nativeは現在の85459項だけから始める。

値はcontent*x^shift*integerPolynomialで保持し、負指数を正確に扱う。
係数と実指数を比較して全new+factor*normalizedPivot=oldを検証した。
parameter座標の単項式だけを反転し、未知normal座標、数値parameter、非単項式の
field pivotは使わない。gcd、Gram、級数、点・素数・新source jetもない。

新6対の元行順は12,171,9,59,38,65、列順は67,12,41,58,70,57。
全てparameter単項式pivotで、新たに用いたparameter逆元はH,s0,t0。
同じ全経路では以前のq0逆元も残り、全五parameterは引き続き独立な係数変数である。

全14266back-add恒等式を確認した。実最大live144498、選択済みforecast最大144590。
native2.139秒、source gate・コンパイル・native全体6.090秒でhard47秒内。
新RAW出力約5.47MBと入力1.88MBの合計も32MiB未満である。

## exact残余と停止の限定

厳密な残余は200×47、5483非零成分、143477項。最大一成分248項。
実指数の下端は(H,L,s0,q0,t0)について(−1,0,−3,−1,−3)で、負指数を
切り捨てたり級数近似にしたのではない。全source行・列IDとalpha/beta制約を保持する。

最後にeligibleな単元は328個あった。そのcheap-fill shortlist16個だけを全予測し、
16個とも150000項上界を超えたため、どの部分も実行せず止めた。残る312個は
調べていない。従って全単元unsafe、全最小pivot、消去対数36の最大性とはしない。
残余rank47は既知rank83と同じ経路の36対からの帰結だけである。

30段receiptと新6段receiptをRHSへ前向き適用し、将来正当化される残余の解を得た後に
逆順代入できる。これは有限のfactored instructionで、全42normal方向の逆や
尽くしの証明を構成したわけではない。
