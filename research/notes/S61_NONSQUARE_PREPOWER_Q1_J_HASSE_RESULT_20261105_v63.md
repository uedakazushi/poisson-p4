# Q1 の既知 J 付値：一次結果と別算術 d-Horner SOURCE

2026-11-05。受理済み source/code 監査後の一回の実行で、k=0の層が606項で非零となった。
従って一次結果は \(v_J(Q_1)=0\) である。商や他因子の値を計算していない。

一次報告 `S61_NONSQUARE_PREPOWER_Q1_J_HASSE_ENTRY_20261105_v63` の SHA は
`0e25e61b3dda79cfba5d6b51aedabd22b4961dd4509b705a87e78a0e3c59e2fa`。
逆順全項の引き戻し零・再生・606項と symbolic a¹⁸ denominatorの読戻しが通った。
peak27234/29336 slots、integer149391・scalar/root29227 products、RAW1504744 bytes、約0.037秒。
これは primaryであり、別算術の検算をこの noteだけで代用しない。

## 独立 SOURCE

元のcanonical4811整数をd次で \(I_l(a,e)\) に分け、
\[
H=I_{18},\qquad H\leftarrow(4+6e)H+(5a)^{18-l}I_l
\quad(l=17,\ldots,0)
\]
と通常のHorner法で計算すれば、最後のHはprimaryのP0と同じ整数多項式になる。
分母は numeric D5¹⁸と symbolic a¹⁸に分けて文字通り比較する。
この別方法はbinomial係数やR-power cache、nativeの三pass convolutionを使わない。
元係数からの演算終了まで606項を参照しない。

全source input4811とH/temporary各1667、metadata16000を含め32290 slots。
全product上限65000、stage-awareL1上限361 bit（preop保守364、hard518）、RAW4 MB。
このSOURCE programはHornerの値を計算していない。
zeroordersはk0のため空であり、606全係数と分母・非零を後続peerで要求する。

既知JがQ1の因子でないという結論は、他の普通因子の互いに素やunit idealを意味しない。
FULLUの十分開の補集合、実bad Fitting scheme、全57/239は未解決のままである。
