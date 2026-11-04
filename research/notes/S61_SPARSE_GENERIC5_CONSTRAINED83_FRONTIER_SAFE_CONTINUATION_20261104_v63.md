# 36段の残余とFULL-safe frontierに束縛した一回の単元継続

2026-11-04。凍結済み36段の厳密残余200×47・143477項を使う。
source-only frontierは旧失敗16位置を除き30候補を完了し、そのcursor29で
元(55,78)=(16/405)t0をFULL safeと証明した。whole coefficient/scratch予測148418、
complete RAW予測13477441bytesで、親の条件付き許可の前提を満たす。
gate report、全row-update forecast、literal候補、同じ実残余のhashを束縛した。

その前提が閉じた後、一回だけ新しい算術workerを起動した。既存36段も別経路も
再実行しない。最初の候補を実残余で再確認・再予測してから実行し、以後は同じ
最大16cheap-fill shortlistだけを完全予測する。global minimumとはしない。

exact Laurent表現content*x^shift*integerPolynomialを保持し、parameter単項式だけを
反転する。各whole pivot前のL+2k+3b+64、support-set scratch、符号付き指数表現、
numeric分母LCM/primitive-content/実exponentによるRAW上界を守る。RAW入力ファイルも
32MiBに算入し、全receipt・残余・次の未完了選択のreserveまで含む。gcd、非単項式pivot、
Gram、級数、点・素数やsource jetは使わない。

新しく完了した単元対は(55,78)の1個だけで、同じreordered経路の総数は37個。
全3534個のnew+factor*normalizedPivot=oldをQQ Laurent演算で直接確認した。
新しい逆元はt0だけ。実最大live148078、選択した全scratch予測148418、
native0.537秒、入力・コンパイル・native全体3.784秒でhard47秒・2GB内だった。
RAW入力3.21MBと新出力約4.67MBの合計も32MiB未満である。

厳密な新残余は199×46、5358非零成分、147370項。元source行・列IDと全alpha/beta
制約を保持し、旧結果も不変である。残余rank46は既知rank83と同じ経路の37単元対
からの帰結に限り、新しいrank試験ではない。

新状態ではeligible単元310個のうちcheap-fill shortlist16個を予測し、全て現在の
150000項forecastを超えたため、どれも実行せず停止した。残り294個は未調査である。
36段の状態で得た他候補のunsafe予測を、この変わった37段の状態へ移してはいけない。
全候補unsafe、37が最適消去対数、または新しいrepresentationが不可能とは主張しない。

36段までのreceiptと本1段をRHSへ前向き適用し、将来正当化される残余解の後に逆順
代入できる。これはfinite factored instructionのみで、全42方向の逆、normal ideal所属、
飽和や一般fibre尽くしはまだ構成していない。
