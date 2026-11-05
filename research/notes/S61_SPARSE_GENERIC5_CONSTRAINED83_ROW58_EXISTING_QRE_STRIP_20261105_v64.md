# row58 の既知 Q/R/E 因子：prepared integer worker

このファイルは固定951項だけの実行プログラムの数学的仕様であり、SOURCE 凍結時点では除算は未実行である。SOURCE 予算は同日 ROW58_QRE_STRIP_PREFLIGHT を参照する。

入力を固定 Laurent 単項式 H^(-1)t0^(-1)、QQ content 1/353894400000、primitive 整数 P に分ける。Q、R、E の順で整数擬除算 a^d P=FV+R を行い、全係数を掛け戻す。試行商と剰余は無圧縮の canonical binary、各試行の入力/出力は content-addressed index として保存する。零剰余時のみ Gauss の補題による V/a^d の整数primitive性を検査して次へ進む。各因子で最初の非零剰余をそのまま保存して止める。最後に得た因子をすべて掛け戻し、初期正規化整数の value hash に一致することを検査する。

F9・未知座標・新しい因子の除算、一般GCD、行選択、行列更新を行わない。出力の指数は既存chart単元の評価を整理するだけである。既存の全168×15と全RHS、他14列、P290一回、すべての補集合はそのまま残る。停止時は完了試行とpartial/failure bytesを保存し、自動再試行しない。数値owned68784、fresh RAW予算、60秒、2GBのいずれかを越えれば停止する。

ROOTレビューSHAを引数で固定した --execute 以外は除算を開始しない。実行後の数学的結論は独立readback/因子検査と実結果reportの範囲に限る。
