# 同じG115/target19をliteral pair行で表すsource-only費用

2026-11-04。新scalar値・minor・rankの計算は行わない。
前の4＋4＋3 adapterと同じ残余一行G115/current11・original32 target19を扱う。

既証のpair-only四行の選択pure六列はrank一である。Bhat4はその三つの
source行の既存polynomial結合であるから、P4の下へ追加するBhat4をliteral
G20/current2のpure行へ替えると、二つのborder determinantの比は変わらない。
これは行の比例とdeterminantの多重線形性だけであり、全17pure列の比例を
主張しない。G20のanchor22は保存済みの非零polynomialで、Saux非零と合わせて
対応D5も一般点で非零となる。D5を実展開して非零検査したわけではない。

g=N[6,1]、h=N[10,1]、rho_i=c_i aFC²として

    Ksmall=g G115/c10−h G85/c6
    D5=det(P4,G20 ; columns22,11,13,25,15)
    D6=det(P4,G20,Ksmall ; columns22,11,13,25,15,19)

と置く。既存B3 Cramer恒等式からKsmallのmiddle-heavy三列は零であり、
D(g D5)上の実残余成分はc10 D6/(g D5)となる。残るpair-heavy消去が加えるpure行は
選択六列で同じ一方向にあるため、このborder比へ寄与しない。
g=0およびD5=0の補集合を保持する。gとD5はFULLUの既知単元ではない。
四つのliteral pair D5 chartは、Bhat4の行結合からD(Saux)を覆う。

計算候補の分子は同じ一成分D6だけに固定する。既存の全108source値、
24pair値、Q9九値をatomとして再利用し、新しいG115/G85 source DAGだけを
必要部分へ制限する。QQ分母、a-Laurent offset、既知B poleを全て明示する。
P4の4minor十五個とbottom2minorのLaplace和を、一つの整数格子へ最後に融合する
有限recipeの支持・total degree・高さ・最後の使用時点を計算する。
これは係数を展開した結果ではなく、実workerの安全なmajorantである。
全七remaining tailとactual full57 source liftを保持し、残余rank四・全FULLU57は未決。
