# μ六 marked H+quartic 九候補の全周囲成分性の独立確認

2026-11-07。PRIMARY E。producer と異なる Fraction 実装で、九つの成分性を独立確認した。
全marked source の分類・全磁気の一軌道還元・既111及び相互相違は、別作者の幾何監査を要する。
本receiptだけで成分数を増やさない。二つのnonclosing source及びD5は再計算していない。

実装は既存 `classification_first_essential_quintic_review.py` に isolated mode
`--mu6-marked-ambient` を追記した。旧十FunctionDefのASTは不変。
旧公開89766547版の全bytesと当時のcurrent jobは
`inputs/first_quintic_reviewer_89766547_frozen_20261107.json` に凍結し、旧八源packetの再構成を保った。
producer・legacy・CASはimportせず、既存の独立Fraction算術だけを使用する。

実行前契約 `jobs/plane_hyperplane_quartic_mu6_ambient_review.json` のgateはPASS。
独立report `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU6_AMBIENT_REVIEW_20261107.json` は
11.157秒、全九件PASSである。
各canonical pointを全六磁気sourceの厳密代入から戻し、全350Jacobiの零、
全350×150の線形化D、Euler25列/GL25列/scalar/磁気六列の実57方向、
Dとの全積零、各非零minorと完全rank、全126実Pfaffian quinticとessential rank三を確認した。
全件でEuler rank24、実作用及び磁気込みsource rank48である。

五件はrankD102、raw tangent48と実軌道48が一致するのでprojective23のsmooth unique component。
残四件はrankD101である。実作用から取った48列Fと、切断座標48行で零となるnormal Nについて、
D N=0、Fの非零48minorを確認した。Dの独立101座標列をe_jとすると、
補正したIFT列はU_j=e_j-F(F|rows)^(-1)e_j|rowsである。
[F,N,U]のdetは列操作で[F,N,e]と等しく、残49行のdetに縮小して
保存された全150直和detと厳密一致することを確認した。
従ってchosen sliceの接線を全周囲と誤認せず、真の102次元cutを使っている。

全350二次ベクトルQ=Jac(N)を独立に戻し、[D|Q]の非零102minorと、
Dの全150列を消す余核functional及び非零pairingを確認した。
IFTで101変数を消した唯一normalの残式は非零二次leadingを持ち、被約cutは原点。
全raw局所次元は48以下、実48軌道があるので等号となり、そのprojective23軌道閉包は全成分である。
smoothnessはD102の五件だけに主張する。高次jetやnil指数は計算せず、分類に必要な最弱支持で止めた。

この全周囲証明と、全marked pair/source及び旧成分比較の独立証明がそろった場合のみ、
対応する九つの相異なる成分を台帳へ採録できる。全G/E及び成分数上界は依然未証明である。
