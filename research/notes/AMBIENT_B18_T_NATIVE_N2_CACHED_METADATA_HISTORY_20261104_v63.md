# N2 cache 計算の progress hash 失敗履歴

2026-11-04。実N2計算後の最終metadata検査で、report保存がphase.jsonを
更新し、その直前に採録したphaseのhashが古くなる問題が見つかった。
全数学恒等式は既に通過しており、算術は再実行していない。
最終版はmutable progressをimmutable evidenceから外した。

元の失敗report・実行source・phase・原因記録の四ファイルを、そのまま
metadata-failed-attempt.tar.gzへ保存した。archive memberごとにhashを束縛する。
内部の失敗reportは履歴であり、現行検証結果として再解釈しない。
最終版のreportと全実係数ファイルは変更していない。
