# K7W 後の標数零の有限 frontier は厳密に二次式

2026-11-04。一次計算の厳密な QQ gcd 証明書。独立レビューは別 gate。

保存済みの monic H15 と、新しい K7W 小行列 N202 の remainder R14 だけを
入力に使った。source、NF、小行列 determinant、点、prime を再計算していない。
入力 report は `S61_NONSQUARE_P100_K7W_FINITE15_RESCUE_20261104_v63`、
SHA-256 `6c6bf67338f63d991114572847886be2118003ce247739bc6beb0e89943ab7ae`。

## 厳密な結果

\[
\gcd_{\mathbb Q[z]}(H,R)
 =g(z)=z^2+\frac{200}{13}z+\frac{22800}{13}.
\]

これは mod32003 の次数上限からの推定ではない。保存証明書は

\[
H=c_H H_{\mathrm{prim}},\qquad
R=c_R R_{\mathrm{prim}},\qquad
g=S H_{\mathrm{prim}}+T R_{\mathrm{prim}},\qquad
H=g Q_H,\quad R=g Q_R
\]

の全 QQ 係数を含み、それぞれを文字どおり検証した。
\(H_{\mathrm{prim}},R_{\mathrm{prim}}\) は positive-leading primitive ZZ
多項式、\(c_H,c_R\) は非零 QQ scalar。従って divisibility と Bézout
の両方向から g が厳密な monic gcd である。

既存の五個の actual 9-minor の正規化後の ideal が `(H)` であり、
新たな N が `HQ+R` なので、この六個の固定小行列の ideal は、許容 P100
曲線の unit ring で `(g)` となる。標数零の固定 cover frontier は
厳密な二次有限 QQ algebra である。H が既知の曲線 unit 因子を持たない
という親 gate により、g にそのような因子はない。

これは実際の quartic 行列の rank failure を意味しない。全 P100 の
pointwise 下限149は維持され、`D(g)` 上では六個の小行列から下限150を
得る。全曲線での150・全239・global 9-target right inverse は未証明。
因子分解や別の source/target による追跡は実行していない。

## サイズと一回の計算

H15 は16係数、R14 は15係数。raw R の最大 numerator/denominator は
8585/7358 bits、primitive ZZ height は H114 bits、R7940 bits。
29次 Sylvester 型の Hadamard bound120754 bits を記録してから実行した。

FLINT `fmpq_poly_xgcd` は一回のみ。native0.219秒、完了 wrapper0.947秒。
最大実測係数14976 bits、peak220 polynomial slots（31 input slots込み）、
serialized117 slots、全 packet約0.73 MB。

最初のコンパイルは standalone C++ の `fmpq.h` / `fmpz.h` の include 不足で
0.413秒で失敗し、native 計算には到達しなかった。失敗時のソースと receipt
を `failed-compile-source.cpp`, `compile-frontier.json` に保存し、二つの header
を追加して同じ一回の数学ジョブを実行した。gcd の再計算はない。

コード:
`code/s61_nonsquare_P100_K7W_finite15_QQ_gcd_20261104_v63.py` と `.cpp`。
出力:
`reports/S61_NONSQUARE_P100_K7W_FINITE15_QQ_GCD_20261104_v63/`。

- report SHA-256:
  `2bcbfb40616f46c7ef98675bd2163d510cfed676c63ea42312ad8ca844272e86`
- native SHA-256:
  `16d69e0a0d38d1696411195ab91cecdf8aa2485093179539f52bda29cbd29e0a`
- source preflight SHA-256:
  `c84545b1906303c9f969d419f9a9fa1cb7328841dab2c94f996dffd8dd950f1d`

全入力と出力を固定する。新しい arithmetic job の自動実行はしない。
