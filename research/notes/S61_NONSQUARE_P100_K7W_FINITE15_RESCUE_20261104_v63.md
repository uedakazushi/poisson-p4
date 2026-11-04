# 新しい K7W 小行列による有限例外の縮小

2026-11-04。厳密な一次計算 receipt。独立レビューは別 gate とする。
全 P100 曲線上の quartic pointwise rank 下限149および pure90 constructor は
維持される。本記録は全曲線の下限150や実際の rank failure を主張しない。

## 固定された一個の小行列

別の source-only receipt に固定した八行 `[9,2,3,5,7,6,1,4]` に、
原319行列の行226、literal normalized C68 times physical \(v_3=W\) を
追加した。九列は old local98 の `[0,1,2,3,4,5,6,9,28]`。
既知の五つの曲線 unit による行 LCD だけを用い、全81係数の clearing
identity と全204個の Bareiss exact division を検証した。旧 K2W 小行列の
計算・データは再利用せず、共通の汎用コードだけを使った。

新しい determinant \(N\) は次数202。既存の monic \(H\), 次数15に対し

\[
N=H Q+R,\qquad \deg Q=187,\quad\deg R=14
\]

を文字どおりの QQ 多項式恒等式として検証した。R は非零で、15 terms。
全15639個の QQ 分母検査が32003-unit gate を通った。N や R の leading
coefficient が unit であるとの仮定は使っていない。

## 有限モデルの結果と厳密な範囲

\(V=\mathbb Z_{(32003)}\), \(S=V[z]/(H)\) は V 上自由15である。
全特殊ファイバー上の厳密な modular gcd と Bézout 恒等式は

\[
\gcd(\bar H,\bar R)=z^2+14786z+21448\quad\text{in }\mathbb F_{32003}[z]
\]

を与えた。従って \(S/(R)\) の特殊ファイバーの長さは2であり、unit gate
は通らなかった。N/R の次数低下や localizer の点を取り除いていない。

それでも、有限 V-module \(T=S/(R)\) の特殊ファイバーは二元生成なので、
Nakayama により T 自身が V 上二元生成である。従って

\[
\dim_{\mathbb Q} T\otimes_V\mathbb Q\leq2.
\]

旧五個の actual 9-minor の共通因子は H なので、この六個の固定小行列で
まだ覆えていない標数零の locus は、scheme length が高々2の有限集合に
含まれる。この有限集合は空でもよく、QQ 上の厳密な gcd の次数2や、
実際の quartic rank が落ちる点の存在を意味しない。
選んだ小行列の非零 open 上では下限150、全曲線では既存の下限149のみ。

既存の mod32003 factor inventory との **literal comparison** では、この G2
は `(z+22576)(z+24213)`、inventory indices1,2 に一致する。これは標数零の
因子分解や geometric point の識別ではない。

## 再現用データと実行範囲

コード:
`code/s61_nonsquare_P100_K7W_finite15_rescue_20261104_v63.py` および `.cpp`。
データ:
`reports/S61_NONSQUARE_P100_K7W_FINITE15_RESCUE_20261104_v63/` の
`source-binding.json`, `native.json`, `certificate.json`, `report.json`。

- report SHA-256:
  `6c6bf67338f63d991114572847886be2118003ce247739bc6beb0e89943ab7ae`
- native SHA-256:
  `a053e50bbe509d9019228ce5ea52fa404ad946c1bf06c765aa96e6fdabd88ba8`
- total4.400秒、native0.122秒。
- peak58545 polynomial terms、うち literal reserve51034。
- serialized1944 terms。出力約1.17 MB。

QQ Bézout、QQ gcd(H,R)、新たな source/column/prime、NF、rank、factor search
は実行していない。失敗した K2W border は再計算していない。
入力・コード・実測出力は固定したまま、独立 peer が別 receipt で監査する。

## 全239列への次の source-only 手順との区別

D0 の独立した SOURCE-only 記録
`S61_NONSQUARE_P100_H15_FULL239_SOURCE_PREFLIGHT_20261104_v63` は、
全33 raw Qhat quartics と full319/H2/H1 の対応を保持する。
one-root finite-degree4 live bound は107743、rank phase bound は81712で、
全33-root 同時保持168679は guard を超えるため採用していない。
この preflight は係数・NF・rank の計算結果ではない。

その report SHA-256 は
`958a02292e9a88ea4151eb3c077f4198b8e5aaf84dcc041fdcabb1cfea691109`。
factor inventory report SHA-256 は
`42e09b5758b7bc7ac0b6bab247b430d498b97d4899a10f1c97923d678b4ff9f5`。
今回の長さ2への縮小は下限150用の固定 minor cover の frontier であり、
full239 用の全 H15 source preflight の対象を自動的には狭めない。
