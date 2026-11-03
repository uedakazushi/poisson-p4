# 第33単位pivotとcofactor定数結合：独立監査

凍結した `S61_BINARY_NONSQUARE_SECOND_SCALAR_PIVOT_20261103_v62` を、**aPFB0CJ≠0、P=ad** 上の指定された一段のrank/cokernel保存操作として受理する。独立検算は2.95秒で終了した。大きいSchur再計算や新しいGBは起動していない。

## 元の二行から実際の59列へ

前の二行packetは、独立監査済みの元26×60行列に対する `M*T=aI2` を持つ。今回、主producerのFraction演算を呼ばず、SymPy QQ[a,e,P]でその元certificateから138項の59-entry transportを再構成した。保存したcofactorの全項と一致した。

さらに実際の25×59 sourceのhashを前packetの出力hashと照合し、metadataが第6行を元第0行へ対応させることを確認した。第6行の全59 polynomialを、certificateの全係数と照合した。

元二行から作る未正規化Schur行

\[
 R_j=a(M_0)_j-(M_0L_7)(M_7)_j\quad(j\ne43)
\]

は、実際の保存行rに対し、全59 entryで

\[
 R=\frac{422053134351}{5039050240}\,a\,r
\]

となることを独立に確認した。この係数は指定開の単位である。従って実際のcofactor像は

\[
 r c=\frac{5039050240}{422053134351}a
\]

で、全多項式等式として再検算した。元の行をentryごとにprimitive化した別matrixは使っていない。

## 五つのcofactorから1を露出する操作

非零QQ定数alpha50,alpha53,alpha54,alpha55,alpha57について

\[
 \sum_j\alpha_jc_j=1
\]

を実際のcofactor多項式で確認した。alpha50≠0である。source列を

\[
 A'_{50}=A_{50}/\alpha_{50},\qquad
 A'_j=A_j-(\alpha_j/\alpha_{50})A_{50}\ (j\ne50)
\]

へ変える操作は、QQ上の可逆行列による列変換である。双対の係数は `c'50=1`、他はc'j=cjとなる。実際の59-entry target行について変換前後の多項式像の一致を独立に掛け直した。

新しい係数1を使ったcolumn50 replacementの行列式も1なので、冗長列を追加する必要はない。そのtarget像は上の単位だから、Schur pivotと指定単位・QQ scalarだけの正規化はA上のcokernelを保存する。

新native sourceの変更箇所も監査した。五つの係数結合の等式を検査してから、列変換前後の全25行imageを別々に計算し、完全に一致することをassertしている。変更の対照は正しく双対であり、行列だけ、またはcofactorだけを変えた操作ではない。

## 別のfresh fixtureと保存出力

新しいnative APIを独立にcompileし、主packetと異なる2×3 fixture

\[
 A=\begin{pmatrix}a&-a&0\\1&e&P\end{pmatrix},\quad
 c=(e,e+1,0)^t,\quad\alpha=(-1,1,0)
\]

で検査した。既知単位はaだけとした。cofactor定数1の露出、column0 replacement、一単位pivotが正常に完了し、残余は厳密に `(e+1,P)` となった。未知因子e+1を勝手に除かないことも確認した。

主実行logの `EXACT_COFACTOR_CONSTANT_COMBINATION 5 EXPOSED_COLUMN 50`、`EXACT_WITNESS_REPLACE 50 TARGET 6`、`TERMS 2755382`、`PASS 1` と、実際の24×58出力SHA256・byte数を照合した。全主input hashも現行ファイルと一致した。大出力そのものは今回再計算せず、保存実行と検査したnative sourceを依存として保持する。

受理するのは計33回のconstructive product-row unit pivot、残余24×58への具体的な還元である。残る八product行・十六genuine cubic行のfull rank、product rank41、full cubic rank57、全cokernelの自由性、quartic充足、ambient除去、J=0への適用は含まない。

独立codeは `code/s61_nonsquare_second_scalar_pivot_review_20261103_v63.py`、reportとfresh fixtureは `reports/S61_NONSQUARE_SECOND_SCALAR_PIVOT_20261103_v63_REVIEW/` に保存した。
