# 元行列から−t0を先に使った別順序の厳密な30単元対

2026-11-04。元236×83、49039項だけから始める新しい一回の算術branchである。
旧24段・29段の結果やunsafe forecastを再実行・変更せず、既存の経路へ足し算した
「累計59段」とも数えない。source-only gate b486f3... に従う。

## 予算と実際の操作

最初は元(192,22)=−t0を明示指定し、その後は全定数／parameter Laurent単項式
候補のwhole coefficient/scratch forecastが最小となる安全候補を選んだ。同点は
元source行・列順。全候補を完全予測するか、単調な部分上界で incumbent に勝てない
と証明してから各段を実行する。最小性の確認が時間内に終わらない段は実行しない。

各whole pivotの前にL+2k+3b+64の係数・作業領域上界を検証した。RAW-byteは
numeric分母LCM、primitive-content桁数不等式、実際の符号付き指数から抑えた。
次の選択が未完了で止まる場合にも全記録を残せるよう、1MiB metadataに加えて
128*remainingRows*remainingCols bytesをreserveした。gzipは使わない。

親の47秒TOTAL上限に合わせ、nativeの次段開始・選択soft guardを39秒にした。
入力保存とコンパイルを含む全wrapperには47秒alarm・2GBを課す。新しい点、素数、
gcd、Gram、非単項式pivot、級数、T/P展開や密逆行列は使わない。

## 完了結果と厳密な停止範囲

別branchの単元対は30個。QQ定数15回・parameter単項式15回である。
全32228個のnew+factor*normalizedPivot=oldを厳密Laurent演算で確認した。
所有係数の最大85851、選択したwhole forecastの最大86226で、150000を下回った。
native39.025秒、入力・コンパイル・native合計41.842秒、RAW約5.74MBで終了した。

残余は正確に206×53、6078非零成分・85459項、最大一成分110項。
152非零行と54零行をそのまま保存する。元行・列ID、正規化pivot行、全target factor、
全行更新forecast、選択証明を各段で記録した。alpha/betaの制約を落としていない。

係数はcontent*x^shift*integerPolynomialであり、負指数も正確な値である。
反転したparameterはH,s0,q0,t0だけ。実指数の範囲は下端(−1,0,−1,−1,−2)、
上端(7,3,6,5,4)。数値や未知normal座標を単元としたわけではない。

最後は39秒の選択guardで止まった。次の605候補のうち125件を調べ、4safe・
120dominated・1time-partialを保存した。この段のpivotは実行していない。
従って「安全候補が無い」「605候補が全て失敗した」とは主張しない。
残余rank53は既知rank83から30単元対を除いた帰結のみで、新しいrank試験でもない。

## 旧経路との差と利用範囲

旧経路は29対で143215項だったが、この異なる順序では30対・85459項となった。
これにより旧項数境界がこの元行列の不可避な境界ではなかったと分かる。
特定の順序の最適性や、残る53方向の消去可能数は未証明である。

任意の236成分RHSへこのbranchの行操作を適用し、将来正当化される残余の解の後に
30pivotを逆順に代入すれば、元source座標へ戻る。これは有限のfactored instructionで、
全42normal方向のleft inverseや残る84coupled方向を構成済みとはしない。
normal ideal所属、飽和、残余の空性・一般fibreの尽くしは未計算。
