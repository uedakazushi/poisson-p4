# H11 d=0 未完計算のCloud引継ぎ

新規研究・再検算を停止。最後の実行 session 93757 は exit 0 で完了。自分が起動した稼働中の計算は残っていない。プロセス全体の ps は実行環境が禁止している。

**保存済み。** H11 c=1,f=−5,L=1,d=0。固定β・固定支持面で第二因子を動かす射影直線を厳密に再構成。旗は既存 S61_MARKED_SIGMA_INTERSECTION_REMOVAL の b3=4/5,b12=b5=1,b6=3,b14=2, 他零。有限Sp4枠で β の非零係数は u4=−(5m²+7)/120,u14=(m+2)/9,u35=1,u49=(5m+1)/9,u50=−5。実Qは比例を除き x(5(m+1)x−6z)∂y∧∂w。m=∞の平方旗を含め全第二因子を覆う。

有限 m の元135の五軌道切断 [5,19,20,23,27] に u12,u13 を残す。QQ[m] Smith分解で定数のみを消去し、(5m+19)(25m−49)/5400 は割らない。自由重みは七つの1と二つの2。F32003[m]で四次561×269行加群が全空間（1325357回、269先頭次数すべて0）。すべての有限 m∈Fpbar で正重み中心環は有限。実際の∞旗の元135中心特殊環ではρ=u5を残し、dim1、ρ0商長76、ρ飽和の根基は標準軸、飽和後ρ0長31、t²=0を F32003 上で再現。Hilbert級数は (1+7T+23T²+35T³−8T⁴−22T⁵−5T⁶)/(1−T)。これで混合標数の周囲消滅はまだ証明されていない。別の torus 極限 x=m,z=1 は正次元中心を持つ別点であり、試験直線の∞には使えない。

**未完の橋。** V=Z_(32003)上で未飽和char0の135式のp-flat closure Yを取る。Sp4の射影座標作用自体がV上で定義されれば、Yはgeneric同変性から安定。元135はYで消え、特殊fiberの計算を上界に使える。閉平方零商DがV-flatならYを経由する。この方針はmps_catalogが概念監査した。具体的なp整性・全相対近傍の証明は未完。

1. S61.reconstruct全係数、vertical_rref、perspective_matrix逆を調べ、V上の群作用を確認。S61_BINARY_PRIME_SPREADINGを再利用。有限chartのSmith変換および群切断minorのp単元性も明示。
2. S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT の明示K・旗接束・Dの平坦性をV上で追う。既存minor (80/21)u5⁴,−4/u5,−u5。Sp4の有限root群または滑らかな局所chartを使う。S61_FLAG_TRUE_CORE_WEIGHTED_REMOVAL のr=−4/25の既存modp重み7,8,9 minor、right inverse、Dの重みを照合。pを含む非正重み底Aに対し実kernelの特殊fiberが有限長ならD-flatと完備NakayamaでA有限。tSはcanonicalな斉次部分加群なので t にある冪零性が従う。指数2を必要としない。
3. 真Dで群補償minor54/25と2、有限Borel4変換、nuisance収縮でρ=1の混合周囲へ運ぶ。実旗端の実加群 t^N S を取り、特殊化と t^N S0 を同一視しない。特殊中心の唯一ρ軸のρ≠0上では前段の混合周囲消滅とGmが使える。実加群の底特殊化が有限長、従って相対A有限、重みで t 冪零。中央Jacobian28とorbit−12/5のp整性を確認。
4. 全P¹_VにYの t≠0 の閉包支持を引き戻し、特殊fiber全点回避とpropernessで標数零の全mへ。支持面条件 P(ker l)=0 は三次式Pの係数に線形な既約4次、一般支持面剰余次数10。l≠0で第二因子はP¹束なので1つの固定(β,Π)の全第二因子が除去できればproper像は16D支持面因子の真閉部分となり、一般β上の全 d=0 を除去できる。現時点でこの結論は認定していない。

**保存場所と再開。** H11_D_ZERO_WIP.tar.gz は work/ の5本の s61_h11_d_zero_*_20261102.py と /tmp/h11-d-zero-pencil-20261102/ の既存入力・途中ログ・終了出力を相対パス保存。ローカルMach-Oだけ除外した。展開後 work/ はrepoへ、h11-d-zero-pencil-20261102 は/tmpへ置く。既実行順 pencil → smith_probe → quartic_smith → code/s61_binary_hyperplane_vertical_rowmodule_20261102.cpp → pencil_flags。M2-binaryを使う。full4がnative用入力、determinant/elimination.logが既出力、actual-infinity-p.m2が実旗端のFp検算。local_samplesは探索用。正式化時はscratchのinspect文字列置換／pickleを独立ビルダに置換する。元135・既存の群作用/旗の依存はrepoのcode,reports,legacyにある。pencil_flagsの実行結果は端末のみなので上に転記した。停止指示以降に検算はしていない。

直近の正式三点組 S61_AFFINE_ANCHOR_RESIDUE_DEGREE_20261102（anchor支持面次数≥10、曲面があればc≥20。T長さ≤19は十分だが未証明）、S61_BINARY_ZERO_PAIR_STRATUM_REMOVAL_REVIEW_20261102 は親へ共有済み。これらと未完のd=0を区別する。
