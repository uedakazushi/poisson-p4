# tail1一個の入力custodyと未実行worker

2026-11-05。`TAIL8_COMPLEMENT_SOURCE`で選んだ列1だけの準備記録。
**T1の列値・κ内積は未実行**。独立source/code監査と親のGOを経て初めて実行する。
八列への自動拡張、別minor・factorization・GCD・点・標数は含まない。

入力は保存済み四κの4920項と、current15四行`[2,4,5,8]`の列1に必要な9268項。
149-node閉包、H2 `[27,31]`、literal anchor22、B一乗のclearingを固定する。
sourceファイルから数値を新たに照会せず、準備済みliteral入力だけをscalar workerへ渡す。
列の全149nodeには和差積しかなく、未知逆元はない。

新しい式は正確に `T1=B*Σκ_i*M[current_i,original32 column1]`。
既に監査済みの九pure列rank一以下とH3零を保つため、同じM12に対して
`A*D15_1=-T1*D12`となる。これは全15 sourceの十分条件へ進むためのsource恒等式であり、
T1自身の非零性、T=0上でT1が非零となる点の存在、実不良Fitting idealとの等号はまだない。

nativeはpacked exponent keyと整数係数の両方を数え、全数値metadata・scalar scratchを
15100 cells確保する。全live term数上限38288により91676 cellsを超えない。
sourceの係数上界359 bits、分母処理395 bits、literal正規化と逆積検査260 bits。
一つのaccumulatorに四積を直接加算し、積多項式や第二accumulatorを保持しない。
検算はforward、逆順で全四積を引いてliteral零、forward再現の三pass。
source積上界1313080対、内積三pass上界42850872対であり、完走時間を保証はしない。

45秒soft/60秒hard、2GB、全RAW 9517864 bytesを保持する。
prepared入力の実サイズを含む全state上界は5852264 bytes。
停止時はJ1KP001一ファイルに全active map、use counts、phaseと演算cursorを保存する。
部分prefixからT1の零・非零を結論せず、自動継続も行わない。
成功出力とprefixの二重全量publicationを予定せず、一つの状態ファイルを使う。
`--inspect-prefix`はstreaming parserであり、再開の承認ではない。

初回準備ではliteral cutとnormalizationを終えた後、nativeの`--help`が
`libflint.so.19`のruntime path欠如で停止した。元code・copy bytes・receiptを保存した。
V2はその四sourceファイルをbyte同一で再利用し、既設FLINTへのrpathだけを補った。
さらにV2初回のmetadata固定値58は誤りで、実際は57 atom recordsだった。
9268+4920という係数項数は変わらない。歴史的ファイル名`literal-atoms58.jsonl`も
変更せず、失敗codeとreceiptを保存して実数57との照合を通した。
いずれの準備失敗でもj1 native係数演算は一度も開始していない。

codeとcustodyの入口は同STEMのV2 Python、nativeは
`s61_nonsquare_tail1_cached_kappa_20261105_v63.cpp`。
値workerのPython入口は`s61_nonsquare_prePOWER_tail1_cached_kappa_entry_20261105_v63.py`。
今回実行したのはcompileと両`--help`のみである。
