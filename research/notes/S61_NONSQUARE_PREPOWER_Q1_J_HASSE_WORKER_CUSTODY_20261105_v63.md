# 共通分母 Q1 の既知 J 付値：SOURCE と未実行 worker custody

2026-11-05。これは独立受理済み canonical4811係数を再利用する準備であり、J の値計算ではない。
\(Q_1=T_1/(a^3C^4)=I_1/D\) の original D197と4811整数 bodyをそのまま用い、headerだけを変更する。
この Q1 は元の S61 の二次変数 q1ではない。

## 普通の多項式環での付値

\[
J=5ad-6e-4=5a(d-d_0),\qquad d_0=(6e+4)/(5a).
\]

係数5aと6e+4の gcdは1（aは後者を割らない）。従ってJは d に関して primitive で、
\(\mathbf Q(a,e)[d]\) の一次既約式から Gauss の補題により \(\mathbf Q[a,e,d]\) でも既約。
これは係数の gcd に関する主張で、両係数が多項式の unit ideal を生成するという主張ではない。
Hasseの最初の非零次数は普通の \(v_J(Q_1)\) に等しい。
5aは既存FULLUの単元なので、未知の分母 chartを追加しない。

固定した \(D5^{18}a^{18}\) を払い、k=0..18について

\[
P_k=\sum_{l\ge k}I_{ijl}{l\choose k}(4+6e)^{l-k}
 5^{18-l+k}a^{i+18-l+k}e^j
\]

を計算する設計である。出力分母は numeric D5¹⁸と symbolic a¹⁸を分けて保存する。
数値の a を選んだり、a=0というFULLU外で評価したりしない。
最初の非零で止め、逆順全項の引き戻し零・再生・全係数読戻しを要求する。

## 全資源と停止

input4811、root cache190、単一二変数accumulator1667。
係数とpacked keyの二セル、全numeric metadata/scratch16000を含む上限は29336。
支持は a≤34、e≤51−k、a+e≤68−kで、各層の上限をSOURCEに列挙する。

入力L1は301 bit。\(5^{18-l}(5+10)^l\le15^{18}\) により全partial/productは372 bit以下、
従来の hard518を維持する。逆順確認も未処理部分和という同じ上限を持つ。
全19次の途中停止位置と二回のrepeatを含めた exact index scheduleの上限は441936 productsで、
先の120万という境界を超えない。45秒 soft、60秒 hard、2 GB、RAW8 MB。

STOP時は入力・root cache・accumulator・denominator/L1・全phase/source/root cursorを
`Q1JHP001` に一度だけ保存し、自動再実行しない。I/O失敗時に二つ目の全prefixを作らない。
準備programはcompileと二つの `--help` だけを実行する。
Jの算術には新しいROOT source/code監査とONE GOが必要である。

本資料は他因子・GCD・tail・商を計算しない。
既存のFULLU十分開 \(D(D_{12})\cap(D(J_0)\cup D(J_1))\) の補集合消滅や、
実bad Fitting idealとの等号も主張しない。
