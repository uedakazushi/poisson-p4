# 全三次残余の二つの元product行：独立監査

`S61_BINARY_NONSQUARE_SCALAR_PRODUCT_ROWS_20261103_v62` を、以下の範囲で受理する。今回の独立検算は1.57秒で終了し、新しいGBや大きいnative Schur出力を再計算していない。

## 元の行との照合と六つの等式

入力は凍結V9の26×60 matrix、SHA256

`6d21ebfbeade3c5a0a7bd2dc091bde5ec274622ddf899add905242164d616ea3`

である。全ファイルhashを照合した後、原matrixをstreamして、第0・7行の全120 entryについてcertificateの全係数と一致することを確認した。entryごとにprimitive化した別の二行matrixを使用していない。

主producerの `Fraction` convolutionを呼ばず、独立したSymPyのQQ[a,e,P] polynomial ringで次を掛け直した。

\[
 M_0 L_0=a,\qquad M_7L_7=a,\qquad
 \begin{pmatrix}M_0\\M_7\end{pmatrix}T=aI_2.
\]

scalar二式と2×2の四式が全て厳密に成立した。L0は64項、L7は34項、Tは206項である。aは指定開 `aPF B0 C J≠0` の単位なので、元の2×60 blockは右逆T/aを持ち、A²へのsplit surjectionを定める。

これは**二行の射影についての分裂**であり、全26行cokernelの零性・自由性・平坦性ではない。product rank41、full cubic rank57、quartic長65、ambient t冪零性も結論しない。

## 最初の列置換と保存native出力

L7の第43成分は非零QQ定数である。第43列だけを全column image of L7で置換する行列Eの行列式はこの定数だから、Eは全source上で可逆である。第7行の該当entryはaなので、続く一つの単位pivotはcokernelを保存する。

native helper `witness-image-stats.cpp` のSHAが凍結V9 manifestに含まれることを確認した。helperの全行image計算、指定単位検査、列置換、Schur計算、known-unitのみの正規化をsourceで確認した。保存logの `EXACT_WITNESS_REPLACE 43 TARGET 7`、`TERMS 1801170`、`PASS 1` と、実際の25×59出力hash・byte数も照合した。今回この大出力自体は再生成しておらず、既受理helperと保存実行を依存として保持する。実行済みnative pivotは計32回という主noteの数え方を維持する。

## 次の行への抽象的transportも独立検算

Tの第0列をtとし、c=(L7)43とする。列置換後の係数w=E⁻¹tは

\[
 w_{43}=t_{43}/c,\quad
 w_j=t_j-(L_7)_j t_{43}/c\quad(j\ne43)
\]

で与えられる。これはQQ多項式ベクトルである。二行それぞれの変換後imageを独立に掛け、aおよび0となることを確認した。さらに、第43列を除く未正規化Schur行

\[
 R_j=a(M_0)_j-(M_0L_7)(M_7)_j\quad(j\ne43)
\]

について

\[
 \sum_{j\ne43}R_jw_j=a^2
\]

を実際の全多項式で検査した。従って第2の単位行は、新しいunit探索無しでtransportできる。これはnativeの前後のQQ/known-unit行尺度をまだ掛け戻していない抽象表示であり、33回目の大native出力を実行したという主張ではない。

## 保存物

独立codeは `code/s61_nonsquare_scalar_product_rows_review_20261103_v63.py`、reportは `reports/S61_NONSQUARE_SCALAR_PRODUCT_ROWS_20261103_v63_REVIEW/report.json`。主certificate、実際の入力、immutable V9、native helper、log、出力のhashを照合した。主packetの全input hashも現行ファイルと一致した。J=0への適用は行っていない。
