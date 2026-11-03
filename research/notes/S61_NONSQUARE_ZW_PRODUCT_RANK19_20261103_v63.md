# 元の Z=W 制限は全域 product rank19：147 を146へ訂正

2026-11-03。ここで扱うのは元の七 light 変数のうち **Z=W=0 だけを課した五変数** である。ell は課していない。残る元の座標は `(v0,v1,v2,v4,v7)`、以下これを `(x0,x1,x2,x3,x4)` と書く。

**全域結果。** 係数開 a*C!=0 で、元の四二次式 `(q0,q1,q3,q5)` の五つの light 一次倍20本は、35次元三次空間中で正確に rank19 を持つ。元の weighted Jacobi chart として用いる開は従来の a*F*B*C!=0 である。d や J を逆にしない。

## 1. 20ではなく19となる恒等式

元の六二次式の全係数から制限を再構成すると、多項式線形式 L があり

`q1=L²/8100`, `q3=L*(eL-400a*x0)/8100`, `L_x1=-40C`。

従って `(eL-400a*x0)q1-L*q3=0` は全 core にわたる非零 linear syzygy である。20-product rank20という仮定は一般点でも成立しない。旧QQ点(1,2,3)の実20×35行列も正確に rank19である。点による全域一般化は行わず、次の二つの全域 gate から下界19を証明した。

## 2. 全 aC-open の coprimality と pair gate

L_x1 が unit なので L を最初の座標とする。q0 は元から x1 を含まないので、L=0 後も `(x0,x2,x3,x4)` の四変数二次式のままであり、その対称係数行列の行列式は

`625a²/1024`。

従って q0 modL は全 a!=0 上で非退化、特に代数閉体上でも既約である。

次に L=x0=0 へ q0,q5 を制限する。x1 は -40C のみで解ける。係数行列の第0行は、単項式順序 `(x4²,x3*x4,x3²,x2*x4,x2*x3,x2²)` で

`(0,(9e+4)/5,5a,-5/2,(21e+16)/10,-5/4)`。

q5 側の clearing と safe row/column 操作は a,C のみを使う。一つの unit pivot 後には5生成元/30項の scalar row が残る。この元の全五多項式に対し、保存された117項・最大次数9の QQ polynomial cofactors が **literal sum=1** を満たす。SymPy QQ 全係数比較と別 FLINT QQ checker で一致を検証し、append+unit pivot の実際の出力も0×5となった。証明・再現に GB、有限体、未知 core 因子の逆元を使わない。

この unit identity は次の小さな消去から構成した。五式のうち三つは `(ad,d,1)` について線形であり、その3×3 determinant が非零定数倍の

`E(e)=19734e²+78099e+50732`

となる。2×2 pivot と d の候補は E と互いに素であるため、E を法として ad,d,a を多項式として解く。四本目の式の remainder も E と互いに素である。この一変数 Euclid の恒等式を元の五式まで戻したものが上の117項 certificate である。

従って q0,q5 の L=x0=0 制限は全 aC-open で独立である。q0 modL が既約なので、q5 modL と共通因子があるなら二式は比例し、さらに L=x0=0 上でも比例してしまう。これは pair gate に矛盾する。よって q0,q5 modL は全域で互いに素である。

## 3. 唯一の linear syzygy

既知 unit による生成元変更で四二次式は `(Q0,Q5,L²,L*x0)` となる。任意の linear syzygy を L で割った商へ移すと、coprimality から Q0,Q5 の一次係数はそれぞれ c0*L,c5*L。L を除いてさらに L=x0=0 へ制限すれば pair gate により c0=c5=0。残りは `b*L+c*x0=0` なので `b=k*x0,c=-k*L` だけである。従って kernel は既知 universal syzygy の一次元であり、20−1=19 が全点で成り立つ。

## 4. 四次の正確な個数：generic146、全点上界146

先行 frozen V49 の actual raw四二次式に対する generic CI 証明を今回も再現した。QQ点(a,e,P)=(1,1,1)で指定四変数線形切断後の次数5 Macaulay56×80に、保存56列の厳密非零 determinant がある。これが nonempty open の regular sequence4を与え、関数体上で barJ=`qfour mod ell` は六変数のCI4となる。全 core の CI は主張しない。

七 light 変数の環 S と J=(qfour)、I=J+ell*(Z,W) を用いる。一般点では barJ がCI4なので J もCI4、ell は S/J の nonzerodivisor。従って

`(I:ell)=J+(Z,W)`

であり、exact sequence

`0 -> S/(J,Z,W)(-1) -> S/I -> S/(J,ell) -> 0`

を得る。次数4で、六変数CI4の次元は48、五変数ZW商の次数3の次元は35−19=16。よって

`HF4(S/I)=48+16=64`、`rank(quadratic×lightquadratic)=210−64=146`

が一般点で成り立つ。係数行列は QQ[a,e,d] 上の多項式行列なので、147-minorは全て恒等零、**全 core の階数上界146** が従う。全点の下界146は追加課題である。H*quadratic の六方向は generic に独立なので、H²を含まない重み4の238方向に対する quadratic products の generic rank は152、全点上界も152となる。

従来の147/153という個数は上の universal linear syzygy を見落としていたため、146/152へ明示的に訂正する。Root の元168×210行列の別QQ点検査でも rank146が確認されたが、global upper bound の証明はこの点検査ではなく generic CI と exact sequence による。

全点CI、全点 quadratic quartic rank146、四次全射、長65、ambient除去は未解決である。重い失敗計算と以前の frozen packet は書き換えない。再現は既存 Cloud activation 後に `python code/s61_nonsquare_ZW_product_rank19_20261103_v63.py`。元全係数、L、全4変数 determinant、literal QQ scalar unit と actual pivots、既存generic CI determinantを検証する。
