# 旧有限残余の16×16多項式行列からのexact eliminant

保存済みk1乗算行列MをQQ定数Dで整数多項式行列Nにする。native FLINTの
polynomial-matrix determinantと、右辺e0だけのfraction-free solveを使う。
全256 adjugateは作らず16個の解Vと非零denominator deltaを保存しNV=delta*e0を
全係数で確認する。deltaと実際のdetの一致は明示して区別する。

既存16本のmonic除算からV(e)=Σvj e^j、Qsum=Σvj qjとして
D*k1*V−D*h1*Qsum=delta。さらにQQ単元のsource normalizationを戻し、元の
h1/k1に対するliteral Bezoutを保存・独立整数演算で確認する。

deltaからQQ contentとaの冪だけを除く。残る非零a多項式は旧有限残余の必要条件。
非定数でも実在を意味せず、未知因子は除去・逆元化・因数分解しない。
残りがQQ単元なら、a非零でh1=k1=0は存在せず旧残余全体を除去できる。

点、新素数、41×41 Sylvester、GB、大きなfinitefiber matrixは使わない。
他の有限Q例外成分や全域rank3、p=0/full57/239の結論は別である。
