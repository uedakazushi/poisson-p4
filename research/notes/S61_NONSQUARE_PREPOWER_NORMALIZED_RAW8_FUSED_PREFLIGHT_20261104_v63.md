# 同一raw8候補の既知column単位除去と融合accumulator前見積

先行raw8 source estimateはlive445860・資料34.34MBで不合格として凍結する。
新scalarは実行しない。このpacketは同じ八行八列への別のexact計算法の前見積である。

四つのpair-only sourceについて、既存全係数証明は
`rho_i=c_i*L, L=a*F*C²`、`Hactual=diag(rho) N diag(1,1,1,1,gamma)`。
heavy column27をL、column31をL*gammaで割ると、各成分はliteral `c_i*N_i,j` になる。
二つともFULLU既知単元なのでp、delta2、R、Sauxを反転しない。raw8detから外す因子は
`L²*gamma`だけである。未確認のpure C因子を仮定していない。

次に各rowのsource exponent下界によるa冪だけを除く。選んだ36pure/8heavy inputと
全atomのQQ/a denominator、全15項のsource cofactor recipeを固定する。
top4 minorとbottom4 minorは同じ15項Laplaceである。

最終の15個のpolynomial積と成長accumulatorのcopyを作らず、一つのdense fmpz係数
bufferへscalar convolutionを直接加える方法を指定する。総degree区間とboxの全格子点を
予め割り付けるので、未知cancellationによってbufferが小さくなるとは仮定しない。
共通QQ分母Dはparameter依存しないinteger。各top/bottom対の分母積がDを割ることを
全15項で検査し、そのinteger multiplierを掛けたscalar積を対応slotへ加える。
各top/bottom polynomialの計算は元DAGのexact演算とlast-useであり、source入力、
scratch、常駐dense出力、15000work reserveを数える。

buffer内の指数は(a,e,d)の有限格子点で、新しいparameter点ではない。原子係数もminor
値も新たに計算しない。scalar convolution回数の上界と、60秒内の完了を保証しない点も
保存する。resource guardが通っても、値の計算には別の明示予算が必要である。
