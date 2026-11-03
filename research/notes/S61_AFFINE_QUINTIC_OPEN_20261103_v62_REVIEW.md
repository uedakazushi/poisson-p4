# 実 monic 五次開と既知五点の norm 分離の独立監査

2026-11-03。公開環境の FLINT 設定を activate した新規 producer と全保存資料は時間以外一致した。
独立の純 Python 係数列で全2056係数の Bézout 等号、二つの商の互いに素性、
非零 resultant16461 と五次の squarefree resultant19651を確認した。
29座標を五次で割った余りへ元の全134式を独立に代入し、全て厳密に零であった。
全五点の元134×45 Jacobian を単項式の直接微分から再構築し、共通 C による
五つの45次 determinant と全 left-inverse、B JZ=0、B Jy=I16 を確認した。

2051正方 high block の kernel は次数4以下の χ,r の組合せであり、五次 gcd で割れるので零。
互いに素な商の syzygy と次数制約から係数も零なので block は可逆である。
一般係数でもその determinant を反転して monic 五次 g を元 ideal 内に作れる。
一般 gcd の次数一定性は仮定しない。全134式と χ の g 余りを保つことで、
実全代数の有限5生成、5×675 presentation が厳密に得られる。

共通斉次次数で Y1≠0 へ移し、16方向を残した有限代数を使う。
y=0 の自由五次環 F0 から実商 C0 への全射は、既監査の普遍既知五点による
generic rank 下界5で kernel の rank が零となる。整域上の自由 F0 の部分加群なので
kernel 自体が零である。有限性や flatness を未証明の大きな商へ仮定していない。
従って全余りは y=0 で係数ごとに零となる。

元 Jacobian の cofactor を直接余りの微分と同一視せず、実 cyclic 同型と
monic g の étale 性を介して JZ dP+Jy dy を余りの微分へ移す。
係数ごとの消滅で余りの s 微分も y=0 では零なので、B の適用により微分は I16となる。
単元分母を払って差商を取れば q=M y、M0=I16である。
Norm(det M) の反転から全 y が残余 ideal に入り、逆包含も成り立つ。
この scheme-level ideal 等号と t norm の単元性により、記載した有限逆行列・norm recipe の
非空開で実アフィン代数は既知の有限自由五点代数に一致する。

一般巨大係数表や具体的な飽和冪指数を計算・保存した結果ではない。
complement の全追加支持、全アフィン五点尽くし、一般次元零は未決として残す。
監査結果: `reports/S61_AFFINE_QUINTIC_OPEN_20261103_v62_REVIEW.json`。
