# nil22 image double_C の全六磁気源を従属三超平面成分へ送る

2026年11月7日。`E_multiple_hyperplane`。以下は全係数計算と明示的な実族による。
別作者レビューを依頼中。新成分は追加しない。

## 実際の多項式族

`NIL22_IMAGE_SOURCES_AND_COMPONENTS.md` §3 の double_C を扱う。
座標は `(q1,q2,x,y)`、全六係数を残した二次式を

    R=m0+m1x+m2y+m3x²+m4xy+m5y²

とする。次を全 ε に対して定義する。

    Q=x(y+1)−1, L=x+1, K=1+εL,
    U=(0,−Q),
    V=(LK, −(1+ε)(y+1)²−(1+2ε)(y+1)−ε),
    A=−(1+ε)y−(2+3ε)−εx.

全非零括弧は

    {q1,q2}=R+Aq1,
    {q1,x}=0,     {q1,y}=−Q,
    {q2,x}=LK,    {q2,y}=V_y,
    {x,y}=0.                                               (1)

`[U,V]=AU`、`det(U,V)=LKQ` が恒等的に成り立つので、(1)は全Rで Poisson である。
全4 strict Jacobiと別の全350 projective Jacobiも直接検算した。
ε=0では元の double_C の全六磁気源へ全括弧が正確に戻る。
族の係数は `QQ[m0,…,m5,ε]` の多項式であり、元の源の係数の反転・特殊化はない。

hで斉次化した実反標準因子は

    Dε=−h(x+h){(1+ε)h+εx}(xy+hx−h²).                     (2)

ε≠0では三超平面は相異なり、その一次式は二次元を張る。
最後の二次式はrank3であり、三超平面の共通平面を含む。
これだけでも既済の従属三超平面＋既約quadricの全支持被覆を適用できる。
さらに以下の全源への実同値により、特定の既存24次元成分へ入ることまで確認する。

## 既済の全八変数源との literal な同定

ε≠0,−1で、実 affine 変換

    X=−ε(x+1),
    Y=y+2+ε(x+1)/(ε+1),
    Z=q2,
    W=−ε(ε+1)q1                                           (3)

を取る。Jacobianの行列式は `−ε²(ε+1)` なので可逆である。
新しい六磁気式を

    Rnew(X,Y)=−ε(ε+1) R(−X/ε−1, Y−2+X/(ε+1))              (4)

と置く。`(X,Y,Z,W)` の順で (1) の全六括弧を移すと

    {X,Z}=f=X(X−1),
    {Y,Z}=Y{(ε+1)Y+2X−1},
    {Y,W}=a(X)Y+f,
    {Z,W}=(ε+1)YW−Rnew,
    a(X)=(ε+1)(X+ε),
    {X,Y}={X,W}=0.                                        (5)

従って `THREE_HYPERPLANES_QUADRIC_DEPENDENT.md` §6,(12) の
**全八変数源そのもの**に、正確に

    a0=ε(ε+1), a1=ε+1,
    (r0,…,r5)=Rnewの(1,X,Y,X²,XY,Y²)係数

を代入したものになる。元の source 実装 `projective_source()` と全六括弧を直接比較し、一致を確認した。
このgeneric開では

    a(0)a(1)=a0(a0+a1)=ε(ε+1)³≠0.

また (4) による六磁気係数の6×6変換の行列式は `ε²(ε+1)⁶` で非零である。
全磁気を保持した全源同値であり、一点又は一部periodの一致ではない。

同ノート第8節と独立reviewが、この全八源の飽和像閉包を
`R4_DEPENDENT_THREE_HYPERPLANES_QUADRIC_24` という24次元全成分と同定している。
従って各 ε≠0,−1 の (1) は同成分に属する。
各固定mについてその閉性を ε=0 へ適用すると、**double_C の全六源と全境界**が同成分へ入る。
ε=−1では (3) のchartを使わないが、元の実族 (1) 自体は正則であり、結論にはgeneric開だけで足りる。

## 証明書と検証状態

producer `code/classification_nil22_double_c_dependent_limit.py` と
`reports/CLASSIFICATION_NIL22_DOUBLE_C_DEPENDENT_LIMIT_20261107.json` は、
全4/350 Jacobi、実因子(2)、全六源端点、変換(3)と逆写像、既済全八源の全括弧、
全磁気変換の非零行列式を確認した（0.645秒）。
有限契約は `jobs/nil22_double_c_dependent_limit.json`。
別作者には producer import に依存せず全式と既済全源を再構成するよう依頼した。

`NIL22_DOUBLE_C_RECOGNITION_FRONTIER.md` のgraded境界・二SS20 pencil・二T02固定衝突の
負結果はこの正しい全源被覆の証明には使わない。高次jet、正式な局所環、未知の平坦性も不要である。
