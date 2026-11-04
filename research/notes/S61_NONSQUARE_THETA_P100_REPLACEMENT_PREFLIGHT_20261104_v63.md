# P100 全曲線判定に向けた既知単元前処理

既に独立検証された P100 の実際の 17×24 行列と、保存済み degree447 の selected
minor だけを使う。曲線のパラメータは

    a=-10(zeta+20)/(9zeta^2), e=(zeta+100)/(9zeta),
    d=-3zeta(zeta+20)/100.

FULLU 上の chart 単元は zeta, zeta+20, zeta+60, 203zeta+2700,
zeta+12 の五つである。これらだけを selected minor から完全に除き、残る
primitive polynomial の exact QQ factorization と全係数の掛け戻しを保存する。
未知因子を単元と宣言しない。

同じ五つの単元について、全 17×24 行列の行・列の共通因子を exact 除算する。
さらに非零 QQ 定数だけで行を整数正規化する。全 408 source entries と row/column
unit factors の積を検証する。A は既存の同じ 17 columns、B は残る七つであり、
新しい selected-minor search をしない。保存 minor の単元による商が expected det
である。新たな determinant、solve、rank、有限体、数値点の計算は行わない。

各行の degree と整数係数 L1 norm から、将来の 119 replacement cofactors の
degree、coefficient-height、ordinary coefficient slots と、最初の候補に属する
17-vector の raw 保存量上限を出す。この raw 見積もりに将来の gcd Bézout 係数を
含めたとは主張しない。上限と source identities が小さくても、全点 rank17 は
まだこの preflight の結果ではない。

将来別途許可された場合、fraction-free solve AX=delta B を実行し、delta が
expected det の非零 QQ 定数倍であることを検証する。残る全 factor を覆う
replacement minor と exact Bézout identity が得られた時に初めて全 P100∩FULLU
での rank17 が結論できる。全 cubic57 への移行は既存の source projection theorem
を用い、全 quartic239 は別問題として保持する。
