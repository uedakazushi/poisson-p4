# 保存1884セルの独立局所恒等式 postpeer

primary `4b1ed027…` を、準備ゲート `9eb3ed13…`、独立worker `492200c9…`、ROOTの明示GO `d2c65bae…` に従って一回だけ検算した。独立報告は `reports/S61_COL53_NEXT17_79_STREAMED_VALUE_20261105_v64_LOCAL_IDENTITY_REVIEW/report.json`、SHA-256 `2b0b2212d91ad5f74481ffc5428a9c33bd1d6cfb018a78c7a3cea7130747dd88`。**全1884局所係数恒等式は PASS**。

元の生産者 evaluator は呼び出していない。独立 bounded Q5 parserで、各保存セルの primitive contentから固定source q_iへの整数係数を復元した。各入力leafは元のliteral operand、既存binary、unit/monomial transform、全5座標へ結び付けた。各COPY/ADD/MULでは保存済みchildrenから別実装・逆順に係数式を計算し、保存済みcurrentと辞書の全係数を厳密に比較した。各段の一致とtopological順序により、全1884ノードに対するsource式の恒等式を得る。

各canonical temporary fileとimmutable stream segmentのバイト列を比較し、全payload SHA、offset、項数、prepared tree digestも一致した。q_iをprimitiveと仮定せず、符号と正整数Gを検査し、空セルでは割り算を行わなかった。source leadと129候補の全5座標を維持したままである。

root2997の必要射影は実際に空で、固定target `(-5,0,9,-3,19)` の係数は零だった。これは H/Lを固定しrelative s≤2とした必要射影の零であり、whole minorの零ではない。新しいpivot、rank、matrix/RHS更新はない。FULL F11を411項のfaceへ置き換えてもいない。

実行1.145秒、数値保持の実ピーク143,670で、事前143,793およびhard150,000以内だった。source整数824bit・一時1921bit・hard2048、RAW27,389,132バイト上界、55秒signal・60秒CPU・2GBを維持した。STOPや再試行はなく、primary・元セル・stream・既存gateは変更していない。新しいsource/IFT/gauge/jet/parameter point/functionalは求めていない。

準備時の構文修正・無関係なparser fixture・メタデータの証明は、prepared gateのhistory/custodyへ保存済みである。本実検算は凍結済みの同じコードとGOを使用した。全1884の証明記録は `all1884-independent-local-coefficient-identities.json` に保存した。
