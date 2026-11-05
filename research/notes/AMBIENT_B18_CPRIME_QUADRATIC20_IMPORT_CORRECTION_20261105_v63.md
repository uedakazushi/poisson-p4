# 4×20 worker の import 修正

V2 は main 呼出前の import で停止した。係数計算・出力ディレクトリ生成は開始していない。元の V2 と SOURCE peer は不変とし、別名 V3 に sparse pair `[index, rational-string]` を辞書へ損失なく変換する局所関数 `rec` を置いた。固定公式・入力・118472 QQ/2048 bit/60 s/2 GiB/32 MiB の上限は不変である。

V3 の `--help` は正常終了した。3 個の局所デコーダ fixture を Fraction と照合した。差分は import、局所デコーダ、当 receipt を実行時 input hash に含める箇所だけである。元コード・修正版の content-addressed snapshot と差分を同梱する。これは計算結果ではなく、SOURCE peer に対する実装由来の補足である。standalone fixture の初回は code/ の import path 指定漏れで停止し、V3 不変のまま harness path だけを補って成功した。
