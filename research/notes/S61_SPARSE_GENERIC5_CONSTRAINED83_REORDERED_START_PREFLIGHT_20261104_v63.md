# 元49039項行列から−t0を先に使う、別順序のsource-only gate

2026-11-04。旧24段、続く5段、全unsafe forecastを変更しない。
元の236×83、49039項行列から始める別の計算経路だけをpreflightする。
新しい多項式係数、pivot、gcd、形式級数は計算しない。

元entry(192,22)は文字通り−t0である。旧経路の最初は(188,37)=2400だった。
元(188,22)=−30L、(192,37)=0であり、新しいpivot行には他の非零成分がある。
従って−t0を先に消去すれば、最初の旧pivot行の他成分は非零な
−A188,22*(row192/(-t0))で変わる。旧literal命令の再実行ではない。
正しい消去の最終row spaceが同じであることは、計算経路が同じという意味ではない。
将来のminimum-forecast経路が何段進むかはまだ計算していない。

最初の一段について、元source支持から全和集合を計算し、行列そのものと
L+2k+3b+64の全coefficient/scratch上界を確認する。符号付きLaurent shiftと
実 exponent を用い、−t0の逆元はparameter係数体の元だけである。
T/Pや旧692953項selector、旧124840/143215項残余をresidentにしない。

RAW上界は数値分母のLCM、primitive-contentの桁数不等式と実指数の文字数から
証明する。全次段で同じguardを使う計画とし、停止時の未完了candidate一覧も
保存できるよう、1MiB metadataとは別に128*(remainingRows)*(remainingCols)bytesを
reserveする。gzipによる上限回避は行わない。

whole workerの安全性は帰納的guardであり、全83方向の完了を先取りしない。
各次の全pivotを実行前に予測し、安全でなければ現在の厳密な残余を保存して止める。
47秒で新しい段を止め、hard60秒・2GB・150000項・RAW32MiB。
定数／単項式全候補のminimum選択が時間内に完了しなければ、その段は実行しない。

このpacketはsource-onlyであり、算術branch、非単項式pivot、gcd、級数、点、素数、
Gramや密逆行列をまだ実行していない。rootによる別の実行許可が必要な段階として
保存する。
