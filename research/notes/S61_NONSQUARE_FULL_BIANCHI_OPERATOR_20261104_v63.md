# Full raw Jacobi/Bianchi operator の小さいsource data

2026-11-04。元135式について、旧full_raw_jacobi_mapsを使う80個・140個の**full** raw Jacobi correspondenceを再検証する。35本のq-leading関係だけをfull Bianchiと同一視しない。新しいquartic matrix、点・field・GB・rank/pivotを計算しない。

元61-variable再構成を現在のbeta coreでspecializeすると、affine bivectorはqのaffine-linear式である。q/t斉次化を `hP=P1(q)+t*P0` とする。q-leadingのquadratic/cubic部分だけでなく、t*L、t*Q_beta²、t*C_beta³も保持する。

80個のy-degree3 Jacobi coefficientsと140個のy-degree4 coefficientsを、旧154 generatorsへ移すconstant mapをfull universal polynomialsで確認する。旧20Bianchi identitiesにより134 retained generatorsへ移し、G6列を零として135列に埋め込む。source correspondence・旧generator indices・全core coefficientsを保存する。

このformal Jacobi trivectorにhPをSchouten bracketし、y-degree4,5,6の35,56,84行を求める。低degree2/3のoperatorは旧normalization/eliminationにより零であることを確認する。全175行について、135個の**実際のq/t斉次化方程式**を代入したcubic polynomialが恒等的に零であることを検査する。degree4行のq-partが旧canonical35と全係数一致することも確認する。

追加の125項関係については、正しいhomogenized identity

`t*G6h - sum K_i(q)*F4h_i = 0`

を別に検査して保存する。tをunitとしてG6を消去しない。

q_fixed=1にdehomogenizeした後、実際のlocal equation weightsとlocal coordinate weightsで各operator rowをweight別に分ける。weight4 piecesから、33genuine quartic columnsのconstant coefficientだけを抜き出した小blockを保存する。これはまだknown-unit pivotや319へのsource削減を行ったものではない。

下流はこのblockだけに一回のbounded known-unit pivotを適用できる。採用した行についてのみfull coefficientsを既知T34/Tn frame transportとimplicit graphに通し、低weight sourceと35unit deletionsへのexact liftを検証する。quarter columnsの独立pivotなら既存quadratic22 kernelと区別できるが、direct Bianchi rowsだけで80方向を得たとは言わない。

指定60秒・2GiB・32MiB output上限を守る。検査が止まった場合は未完了sourceとlogを保存し、同じfull calculationを再起動しない。Math claimは保存されたfull polynomial identitiesとactual row countsに限る。
