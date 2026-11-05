# Quadratic20 worker V3 の限定 import 修正監査

2026-11-05。SOURCE 監査 e7b749ee は V2 の数式と AST を確認したが、import 先が `rec` を export しない点を見落とした。V2 の起動は main・出力 directory 作成・値計算の前に停止した。失敗、V2、SOURCE report をすべて保存し、既存 byte を変更しない。

独立監査は V2/V3 の全文差分を固定し、変更が不足 import の削除、局所的な lossless sparse rational decoder の追加、修正証明書の入力 pin だけであることを確認する。三つの decoder fixture を検査し、V3 の `--help` を独立に起動して import の成功と値出力 directory が作成されないことを確認する。main は実行しない。

数式、source、scope、全114/164検査、118472 QQ・2048bit・60秒・2GB・32MiB 制限は変わらない。root が許可した一回の値計算を V3 で続行する条件を満たす。cubic や追加 source の許可ではない。
