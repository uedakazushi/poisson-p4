# CA原点 self-flag 四次jetの独立検証

2026-11-03。対象は `AMBIENT_CA_ORIGIN_SELF_FLAG_4JET_20261103_v63` の有限四次jetだけである。
producerを再実行せず、独立の有理数・Laurent疎多項式演算と実座標変換で監査する。
後続の有限多項式pair構成はこの監査に含めない。

生のCAテンソルからdivergence補正を再計算してcanonical 150×20基底を復元する。
free10とrelative flag全9方向を照合し、正規化した第一flag方向
`M=E10−2E03` を用いる。18補正方向の150行列を実pushforwardの一次微分から
再構成し、rank15と三次元の全kernelを両側inverseおよび完全行列等式で確認する。
三次で残った三自由度はすべて保持する。

主検証は `g P(g^-1 X) g^T` の直接計算である。群の逆行列を
`I−N+N²−N³+N⁴` で作り、二次monomialの逆座標代入と外積作用を順に計算する。
係数環は `QQ[t,t^-1,y0,y1,y2][epsilon]/(epsilon^5)`。
producerのLie/log作用も保存済み三次intrinsic行列も再利用しない。
三次intrinsicから補正を解き直し、全150係数について次数0–3のmembershipと
四次の消去前・消去後係数を照合する。新しい18補正のimageも同じ完全rank15行列である。

全四次cokernelの定数係数rankは4、独立行は零始まり46,63,72,95。
全150行がこの四行の定数線形結合となる完全係数行列等式も確認する。
独立式を `q1,q2,q3,q4` とすると

    y2=(3t/2)q1−(2t/3)q4,
    (t²−1)²=(2t/3)q2,
    4t²y0+(7t²+3)y1=−2t q3

および各qの逆向きの明示倍数表示により、正確なLaurent idealは

    (y2,(t²−1)²,4t²y0+(7t²+3)y1).

従って先行補正のschemeは `QQ[t,t^-1,y1]/((t²−1)²)` である。
`t²−1` は非零で平方零。CRT idempotent `(2+3t−t³)/4` によって二つの
二重A1へ分解する。各符号 `s=±1` で `t=s+e,e²=0` とすると
`y0=−5y1/2+(3s/2)e*y1`。これに新しい四次補正のkernel三自由度を加えると、
**この正規化された**jet解schemeは二つの二重A4となる。
自動的にfiber内に入る最上位source係数は零に固定するか自明な自由因子として加えられる。
このA4という数を無正規化の全raw jet-moduliの次元と解釈しない。

両符号 `t=±1,y0=y1=y2=0` について補正を独立に解き、群・第一source・第二sourceの
全保存係数を照合する。それぞれ実pushforwardの150×5係数と全次数membershipを直接確認する。
`t=1` の新source free10は `(0,0,0,0,0,0,2,0,0,1/3)`、
`t=−1` では新18補正は零で、いずれも四次flag補正は零である。

これは有限四次の必要十分条件と二つの有限witnessの承認である。
全次数の形式的・代数的pair branch、ordinary以外の全chart、個別CAのconductor、
原点の全odd idealについての結論は含まない。

再現はCloudの既設環境をactivateし、次を新しい出力先に対して実行する。

    python code/ambient_CA_origin_self_flag_4jet_review_20261103_v63.py --output <new-output-directory>

独立コード、入力・証明note・証拠のSHA256、成功状態と時間は
`reports/AMBIENT_CA_ORIGIN_SELF_FLAG_4JET_20261103_v63_REVIEW/report.json` に保存する。
