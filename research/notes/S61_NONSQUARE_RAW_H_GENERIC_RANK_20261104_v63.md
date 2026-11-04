# Raw Hのgeneric rankを7×3だけで厳密に判定

既存14×5のrow2,row4とcolumn0,5のdeterminantは、文字通り
`-(6/5)*B*p²`、p=M[2,0]という非零多項式である。この二行は他三列が零なので、
QQ(a,e,d)上でその二列を消しても他三列は変わらない。他三列のsupport行は
厳密に7,9,10,11,12,13,14だけであり、rank問題はこの7×3へ移る。

この小行列だけに最大三段のBareiss消去を行う。各段の非零pivotが元小行列の
選択minorであること、前段pivotでの割り算がexactであることを確認する。
終了後は選択minor Deltaとadjugateを使って全7行の座標を記録し、全21個の
cleared polynomial恒等式でupperを証明する。残るcolumn kernelも多項式で保存
する。点や有限体を使わない。

元raw H16×7のgeneric rankは、既存の二単元方向を加えた
`2+2+rank(7×3)`である。pやDeltaはgeneric係数体でだけ可逆にする。これを
FULLU全域のrank下界、既知単元、自由source表示へ移してはならない。generic
rankが7未満なら、raw H-linear全7方向を全域で埋める十分条件は不可能になるが、
quarticも含むH1射影の全射性は引き続き別の問題である。
