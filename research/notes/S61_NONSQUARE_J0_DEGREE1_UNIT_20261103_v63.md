# J=0 上の一次多項式係数による第37単位消去

## 証明した範囲

`J=5ad-6e-4=0` とし、以下の環で実際の cubic coker を考える。

\[
 A=\mathbf Q[a,e,1/(aF B_s C(3e+2))],\quad
 F=500a+90e^2+3e-24,\quad C=100a+18e^2-9e-2,
\]
\[
 B_s=300ae+54e^3+9e^2-24e-4,\qquad
 d=\frac{2(3e+2)}{5a},\quad B=\frac{2B_s}{5a}.
\]

前の第36単位消去で保存した **21×55** 行列の第4行（すべて0始まり）について、39個の非零係数、合計45単項式からなる一次多項式係数 \(q_j\in\mathbf Q[a,e]\) を構成し、**元の行の55個の多項式を用いて**

\[
 \sum_{j=0}^{54}q_j M_{4j}=(3e+2)^2
\]

を独立の sparse `Fraction` 演算で証明した。係数と55個の入力多項式は `reports/S61_NONSQUARE_J0_DEGREE1_UNIT_20261103_v63/identity.json.gz` に完全に保存してある。これは関数体上の数値的な rank 試験ではない。右辺は指定した開集合で単位であり、J や未知の因子は逆転していない。

係数の中に非零有理定数があるため、保存した係数で一列を置き換える列変換は A 上可逆である。その単位を pivot として一行一列を除き、同じ cokernel の **20×54** 行列を得る。保存した witness による native 消去を独立に再実行し、全出力が以前に保存した行列と byte 単位で一致することを確認した。最初の57行基底からの安全な単位消去は合計 **37回** になった。

## 探索と再検算

探索では各 entry から既知の単位因子と非零 QQ scalar だけを除き、残りの多項式を `1,a,e` 倍した165列の QQ 係数行列を作った。rank は73であった。既知単位の積を探す過程で `(3e+2)^2` が得られた。除いた因子を戻した後、native code 自身が **元の entry 全体との多項式等式** を検査する。独立 verifier は探索行列やその RREF を信頼せず、元の55多項式と保存係数だけを読んで上の等式を再計算する。

新しい方法が定数係数探索の単なる繰り返しではないことも小例で確認した。行 `(e,ae-1)` と既知単位 `a` について、旧定数係数探索は拒否するが、新しい探索は係数 `(a,-1)` による等式 `a*e-(ae-1)=1` を出す。この小例は探索機構の動作確認であり、実際の21×55行列についての証明は別に上記の全多項式等式で与えている。

## QQ scalar による保存サイズの改善

消去直後の20×54行列は237,849,698 bytesであり、その有理係数の分母・分子の共通部分を行、続いて列ごとに整えた。これは非零 **有理定数** による対角変換のみである。全351,404係数について、別の `Fraction` 実装で

\[
 N_{ij}=r_iM_{ij}c_j,\qquad r_i,c_j\in\mathbf Q^\times
\]

を検算し、整数係数の130,953,506 bytesの行列を保存した。core 多項式を約分した操作ではなく、QQ 上の cokernel は同じである。これらの大きい数値定数が特定の素数で単位だという主張はしていない。

raw SHA256:
`80179a8145d79c1046764da2aabc0c51dfed774e67cdd8d0b969b697732c2c00`

scalar 正規化後 SHA256:
`55174c163211b4252a9f1ccff04de204b5aece4a1c26aab2a539f8eb5d489bb6`

巨大行列は安価に再生成できるので公開 report には重複保存していない。元の行列は前の `S61_NONSQUARE_J0_SECOND_UNIT_20261103_v63` packet の再現手順で得られる。

```sh
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_J0_degree1_identity_20261103_v63.py \
  --replay reports/S61_NONSQUARE_J0_DEGREE1_UNIT_20261103_v63/identity.json.gz
```

元の21×55行列が存在する作業領域では、`--output <新規directory>` により元の行の照合、保存 witness の消去、QQ 対角変換の全係数検算まで再現できる。公開コードは `s61_nonsquare_J0_degree1_unit_span_20261103_v63.cpp`（探索）、`polynomial_one_witness_pivot_20261103_v62.cpp`（実際の一単位消去）、`polynomial_QQ_scalar_primitive_20261103_v63.cpp`（QQ scalar 変換）、および独立 verifier である。

## まだ証明していないこと

残り20×54行列の cokernel が零であること、J=0 上の cubic 欠損集合が有限であること、欠損点での quartic coverage や ambient nilpotence は、この packet の結論ではない。以後の探索の timeout や定数係数探索の拒否を、未知の因子の存在証明と解釈してはならない。Cloud の共有ライブラリ環境を有効にする前に一度起きた loader 失敗は数学的計算に含めず、その後の正常な実行だけを用いた。
