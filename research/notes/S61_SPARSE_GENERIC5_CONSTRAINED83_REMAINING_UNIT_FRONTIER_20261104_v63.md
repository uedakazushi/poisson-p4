# 36段後、既失敗16位置を除く単元frontierのsource-only scan

2026-11-04。入力は凍結された200×47・143477項の残余のみ。既に失敗を保存した
正確な16個の元行・列位置をhashで束縛し、候補列挙の最初の16個と一致することを
確認して除外した。残る312個を同じ決定順で一度ずつ予測する。

row operation、pivot、gcd、級数、多項式係数の積・差や新しいsource jetは実行しない。
既存QP係数のnumeric contentとsupportから、whole-pivot更新の全scratch/RAW上界だけを
評価する。unsafeはその予測方式で実行前に上限へ達する証拠であり、代数的rank失敗や
他のrepresentationの不可能性ではない。

各完了forecastを個別ファイルとして保存し、対応するimmutable cursorも保存した。
native32秒soft・入力保存とコンパイルを含むhard40秒、2GB・150000所有係数・
RAW入力と全出力32MiB。各次記録を保存できる2MiBも先にreserveした。

一回のsource-only scanは残る候補のcursor0〜29を完了した。最初の29個は現在の
whole forecastでrigorous unsafe、cursor29の元(55,78)=(16/405)*t0はFULL safeだった。
その全row updateの予測は148418項、final logical support上界147646項、
complete RAW上界13477441bytes。FULL safeが見つかった時点で即座に止め、
算術pivotは実行していない。次cursor30、未調査は282個である。

既失敗16個と新29個がunsafe、1個がFULL safeという限定であり、全328 unsafeとも
global minimumとも主張しない。選んだparameter単項式と全五変数は同じK5上のもの。
chosen-candidate.txtは元位置、QQ係数、実Laurent exponentをliteralに保存する。

所有する多項式係数数は実際に始終143477で不変。native0.572秒、入力・コンパイル・
scan合計3.979秒。全完成forecastとcursorを保存し、元36段の結果は変更していない。
このFULL-safe gateは親が条件付きで許可した別の一回の算術継続に使えるが、
本packetそのものはsource-onlyである。
