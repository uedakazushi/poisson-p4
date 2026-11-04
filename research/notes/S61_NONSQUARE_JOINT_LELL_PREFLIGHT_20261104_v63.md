# joint (L,ell) frame：monic leading term による損失のない前処理

**四つの二次を正確な単項式 L²,Lx0,ellZ,ellW とみなす案は棄却する。**
最初の二恒等式は Z=W=0 上だけであり、full source では非零 Z/W 補正を持つ。
これらを落とすと既知の二次倍行 rank41 に反する余分な linear syzygy を作る。

old light `(x0,x1,x2,x3,x4,Z,W)` に対し
`r=7*x2/2+3*x4`、new lightを `(L,ell,x0,x3,r,Z,W)` とする。
`ell-2J*L/1215` の x2,x4 係数 D2,D4 は
`3D2-7D4/2=-80ad²/81`。L_x1=-40C なので Jacobian は
`delta=3200*a*C*d²/81`。FULLU の既知単元である。
逆線形写像の分子行列 N について forward*N=N*forward=delta*I を全98係数で
検証する。座標分母は a,d,C だけである。

full実六二次を可逆に

```
p1=8100*q1, p2=8100*(e*q1-q3)/(400a),
p3=ellZ, p4=ellW, p5=q0, p6=q5
```

へ変える。降順 lex で最初の四行は monic leading term を
L²,Lx0,ellZ,ellW と持つ。前二行の全 Z/W 補正を保存する。
同じ重みの積行には26（重み3）と98（重み4）の unit triangular square が
あり、これを正確な source/target Schur pivot として使える。
GB・初期ideal等式は主張しない。

重み3は target91から65へ、four-quadric source28から overlap2を残す。
W*p3-Z*p4 は恒等的零だが、x0*p1-L*p2 は非零で、full係数と分母を保存する。
q0/q5 の14倍行と genuine16を合わせ、正確な31→65入力の circuit recipe を得る。
重み4は target239から141、source319から221へ移る。
残る18 overlap と58 pair倍行、112 cubic倍行、全33 genuine quarticを保持する。
65や141は選んだ unit square の補空間の数であり、実商の Hilbert function
ではない。full239・全点cubic57は未解決。

初段は六二次を多項式分子・既知単元分母へ実際に変換し、係数増大を測る。
G16 と33quarticの新しい parameter 係数展開はせず、既存正確な source と
N/delta による lazy coefficient recipe を保持する。
次の finite-circuit stage で98 pivotと既知38方向を輸送する場合にも、四つの
単項式idealで割ったと誤認せず全 overlap 行を残す必要がある。
