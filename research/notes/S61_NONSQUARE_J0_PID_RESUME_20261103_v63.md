# J=0 PID 継続：保存済み12段を確認、正確な24行70列を保存

2026-11-03。先の
`S61_NONSQUARE_J0_PID_20261103_v63` の「11段完了＋12段目保存済み」を、今回の再実行で
**12段全て確認済み、正確な24行70列行列を export 済み**へ進めた。以前の checkpoint は
その時点の記録として変更しない。

## 確定した進展

Cloud 復旧時に最初に確認したところ、旧 `checkpoint-replay.json` と `resume-matrix.sing`
は存在せず、log は空だった。そこで保存済みの cofactor だけを再適用した。最初の12段の
row-gcd 探索や traced Smith 計算は繰り返していない。

最初の60 CPU秒の prefix replay は final PASS/export に届かなかった。その後、新しい
総600秒・AS 2.2 GBの継続 budget を設け、各段の完了 marker と途中状態の export を追加した。
168.71秒で12段全ての exact replay と24行70列行列の export が完了した。log の各段と
最終 `EXACT_SAVED_12_STEP_PREFIX_REPLAY_PASS`、`SAVED_STATE_12_EXPORT_COMPLETE` を
実際に確認した。arithmetic error はない。

確認した全12段は、\(\mathbb Q(e)[a]\) 上での cokernel 同値である。polynomial variable
\(a\) について除く因子は既知の \(a,F,B_s,C,3e+2\) のみであり、係数体に由来する非零
\(e\)-多項式の除算と区別する。\(J\) は反転していない。

保存した正確な行列は

```
work/nonsquare-J0-pid-v63/continuation600-v1/replay-1/state-after-12.sing
```

で、24行70列、43,952,280 bytes。`ring R=(0,e),(a),dp` と既知単元の定義を含む
self-contained な Singular 入力である。最初の12段を再探索せず、この行列から直接再開できる。

## その後の実際の状態

export 後の残り約431秒を使い、この24行から row-gcd を継続した。新しい一つの pivot の
cofactor は保存され、log には `EXACT_ROW_GCD_PIVOT_e 1 LEFT 24 ...` が残った。
しかし、その表示は最終 matrix 更新より前であるため、**新しい13段目の matrix 更新の
完了までは確認していない**。

その後 Cloud execution environment が再び利用不可となった。復旧直後の確認では
Python/Singular process はなく、600秒 controller の `report.json` と継続アルゴリズムの
`receipt.json` は両方とも存在しなかった。従ってこれは「600秒 timeout が記録された」と
読み替えず、**最終 receipt 前の環境中断**として保存する。新しい一つの witness は未確認の
途中データとして残し、確認済み24行行列に混ぜない。

この継続でも、全余核の消滅、実際の欠損曲線、J=0 の full cubic rank-defect support の
有限性は証明していない。既存の product-rank41 欠損有限性を full cubic の結果に置き換えない。

## 公開保存と再開

`reports/S61_NONSQUARE_J0_PID_RESUME_20261103_v63/confirmed-state.tar.gz` は、確認済み
24行70列行列、その replay script・log・receipt、未確認の新 cofactor と継続 log を元の
`work/...` 相対 path で保存する。`report.json` は全 member hash、確認済み／未確認の境界、
中断理由を固定する。実行 binary は含めない。

元の36行70列入力、全2520成分の QQ 係数検査、13個の bivariate known-unit pivot、二つの
Smith timeout、backend failure の除外は、前の公開 checkpoint に保存済みである。
今回の archive はその確認済み12段以降の delta である。

次の操作は archive を repository root に展開し、まず確認済み `state-after-12.sing` から
始めることである。未確認の新 witness を使うなら、その pivot 行の Bezout 恒等式と既知単元
性、および matrix 更新を先に再検証する。環境中断前に計算済みと仮定しない。

現状では rational-function coefficient の matrix 更新が主な費用になっている。保存された
行列から係数の \(e\)-denominator を行ごとに払い、native bivariate polynomial arithmetic で
既知単元・係数体単元だけを正確に取り除く方法も、変更のない Smith 再実行より検討価値がある。
この方法はまだ実装・証明した結果ではない。全 J=0 の有限性と quartic/ambient の除去は
引き続き未完了である。
