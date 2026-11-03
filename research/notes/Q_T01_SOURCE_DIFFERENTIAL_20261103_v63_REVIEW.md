# 全 Qa における actual T01 特異性の独立監査

2026-11-03。新規 producer report は時間と gzip 圧縮時刻に伴う hash 以外一致し、
圧縮を解いた全行列 JSON は完全一致した。圧縮済み公表資料自体の hash も照合した。
閉じた Fréchet 和の helper を使わず、双対数 (N+εE)^k、ε²=0 の反復乗算で
全25行列方向を微分し、六係数方向と合わせた実150×31行列を独立に再構築した。
両表示の厳密 QQ(s) 行列、各 rank24、和 rank30、小行列式を再確認した。
直接単項式係数の GL5 作用は rank17で、両表示の equivariance 恒等式と一致した。

5個の generator 変更では、最初の25成分 N^r を固定し、六係数方向の補償を
150式から独立に解き、保存 kernel と一致した。strict stabilizer の8方向、
conjugation の3方向との和の rank7、regular nilpotent centralizer の5方向も確認した。
polynomial generator change の可逆4×4変換と滑らかな標数零 stabilizer の作用は
実際の fiber 内に積分でき、合わせた微分 rank7 は fiber 次元7以上を与える。
fiber の embedding dimension も31−24=7なので、一般 source fiber は正則7次元である。
target の正則性や特殊 a の source fiber の正則性をここから結論しない。

共通18次元空間は実 GL 軌道17方向と Qa 曲線方向に一致する。
両表示は全域で actual reduced T01 に写るという別の独立監査済み結果があるため、
二つの微分像の和30次元は reduced T01 自体の接空間に含まれる。
T01 の affine dimension26より大きく、一般 Qa は特異である。
接空間次元30以上の集合は閉、a↦Qa は全 A1 上の非零 tensor の射なので、
下界30は a=0 と cover 分岐値を含む全有限 a に成立する。
非零 cone 点の自由 scalar 因子を除くと projective tangent の下界は29となる。

二つの source 値や像の接空間の違いから target の局所枝の本数や
normalization fiber を推論しない。全 T01 completed ideal も決定していない。
監査結果: `reports/Q_T01_SOURCE_DIFFERENTIAL_20261103_v63_REVIEW.json`。
