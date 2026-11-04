# P100 low98：cached input と NF last-use の細部

新しい係数値は計算せず、受理済み十個の pure cubic kernel と六個の pivot cubic の full physical84係数を読む。全16×84を合わせた保存 RF support は16076 terms。cached numerator degree の最大は13、denominator degree は11。84列ごとの degree/support、各行の値の有無を保存した。これは元16cubicの raw tails の全てが既に同じ frame で実体化されているという意味ではなく、kernel＋pivotという source basis の入力である。

受理済み pure98 NF linear DAG の needed nodes を固定した順序で last-use 解放する schedule も保存した。一次 value の peak は121。この数は一つの RF value が何 termsへ増えるかとは独立である。scalar multiplier と joint-frame transport の RF degree growthを含む bound が必要なので、121だけを150000 coefficient supportの証明へ読み替えない。

この細部は単純な全70行保持454974という既存 forecast を、streaming型の実装で改善する入力である。ただし全98 polynomial matrixの degree、Bareiss state、global determinant、gcdは依然として未計算・未認証。coefficient/NF/rank/点/primeの新しい job は実行していない。
