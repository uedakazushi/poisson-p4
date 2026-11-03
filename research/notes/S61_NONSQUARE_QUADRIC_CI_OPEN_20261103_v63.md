# 次数9の明示的な開で actual 四二次式は完全交叉

2026-11-03。raw六 light 変数の四二次式を `Q0,Q1,Q2prime,Q3` とし、既知開を a*P*J!=0 とする。先行 actual source と binary reduction の全恒等式を再検証し、**追加の明示的 chart 条件 R_vertex!=0 上でこの四式がCI4** であることを証明する。R_vertex は自然 core `(t,v,H)` における次数9の多項式で、未知因子のまま保存する。R_vertex=0 の場合は解決していない。

## 1. 四変数の determinantal presentation

既存の厳密表示は `c=6075` として

`Q1=AB/c+W*C1`, `Q2prime=40X*B/c+W*C2`, `D=40X*C1-A*C2`。

六変数の可逆座標 `(A,X,W,B,S,U)` の inverse determinant は `1/(200aP³J)`。C1,C2 は最初の四座標だけを含む。三式 Q1,Q2prime,D/c は `[[A/c,40X/c,W],[C1,C2,-B]]` の符号を変えた2-minorである。

また B=W=0 上の D は binary quadratic で、その X²係数は `1280a³J²/243` という既知 unit である。

## 2. 基底は純一次元の projective curve

Q1,Q2prime は共通因子を持たない。共通 quadratic 因子なら比例するが、W=0 では AB と XB であり A,X は独立なので不可能。共通 linear factor F が W に比例することも不可能。W=0 上で gcd(AB,XB)=B なので、F は適切に定数倍すれば `B+kW` となる。一方 `W*D=40X*Q1-A*Q2prime` だから、F と W が互いに素であることより F は D を割る。これは B=W=0 上の非零 binary D に矛盾する。

従って `(Q1,Q2prime)` は height2。三つの2-minor ideal も height2であり、2×3 maximal-minor ideal の expected-height theorem / Hilbert–Burch により、四変数 `(A,X,W,B)` の quotient はCMで、projective support C は純一次元の次数3 curveである。C の既約性やreducednessは仮定しない。六変数での determinantal locus は C の vertex `V=P1={A=X=W=B=0}` を持つ cone である。

点ごとの恒等式から

`V(Q1,Q2prime)=V(B,W) union V(Q1,Q2prime,D)`

が成り立つ。W!=0 では D の定義から従い、W=0 では B=0 または A=X=0 である。

## 3. Vertex binary resultant と R_vertex

V上でQ0,Q3を制限した二つの binary quadratics は、先行 actual pair chart で B=0 を課したもの、すなわち係数列 `[2,4,5]`（light monomial vlight²,u*vlight,u²）である。

自然 core `t=a/P, v=3e+20t, H=P+8t` と既知 unit row scales/row subtractionを用いる。固定4×4 Sylvester determinant の厳密恒等式は

`Res_vertex(Q0,Q3) = -729*t*R_vertex(t,v,H)/4096000000000000000`。

R_vertex は零でない次数9の多項式で、全係数を certificate.json に保存した。原の F,B,C の積であるとは主張しない。t は既知 unit なので、R_vertex!=0 と「vertex上でQ0,Q3に共通零点がない」は同値である。source replay は元804項からこの binary pair を再構成し、固定次数の全Sylvester determinant を照合する。

## 4. Proper projection による height4

R_vertex!=0 とする。四二次式の共通零点を Y⊂P5 とすると Y は vertex V と交わらない。

第一の部分 `Y∩V(B,W)` は P3 内の閉集合であり、その P1⊂P3 の vertex と交わらない。projective dimension theorem により、この部分の次元は1以下である。この議論により、この開では別の gcd-exception の除外を追加する必要はない。

第二の部分は cubic cone 上にある。Vを中心とする射影を Y に制限すると Cへの regular proper mapとなる。Cのどの幾何点上の fiberも、その ruling P2の閉部分であり、vertex line P1と交わらない。projective dimension theorem によって各fiberは零次元以下である。したがって射影はproper quasi-finite、特にdim(Y on cone)≤dimC=1である。

両部分を合わせてdimY≤1。六変数多項式環中の四式のheightは少なくとも4であり、四生成元なので正確に4。従って regular sequence4となる。有限体、標本点、base C のsmoothnessを用いた推論ではない。

## 5. 元の六二次式の四次階数への帰結

`a*P*F*B*C*J*R_vertex!=0` では ell chart が可逆で、上の raw CI は元の四二次式 mod ell のCIと同値である。元の五変数Z=W制限の20-products rank19は別の全域QQ packetで証明された。従って `(I6:ell)=J+(Z,W)` と exact sequence から、light次数4の商は `48+16=64`、quadratic-products のrankは正確に146となる。H倍の六独立方向を加えた weighted quadratic-product rank は152である。

これは二次式の四次倍の階数である。全三次の点別rank57、残る quartic quotient の全射、横断長65やambient除去は、このCI-openだけからは結論しない。R_vertex=0 は正当な未解決のcase equationであり、全 core の既知 unit として扱わない。

再現は既存 Cloud activation 後 `python code/s61_nonsquare_quadric_CI_open_20261103_v63.py`。先行二つの actual source replay、binary D のunit係数、四座標frame、固定4×4resultantを検証する。GB、追加finitefield、大きな四次行列展開は行わない。
