# nil22 quadratic_Aの全七変数源は実symplectic曲線でMPS40へ入る

2026年11月7日、E_multiple_hyperplane。producerによる全源の証明、独立レビュー待ち。
**全六磁気係数とλを持つquadratic_A源は既知MPS40の閉和に含まれる。**
新成分性、停止した十八変数局所idealの支持、全形式変形の分類はこの証明に不要である。

## 1. 全六磁気と全λを保持する実座標変更

α=λ−1/2、β=−λ−1/2、I=m0+m4、
R=m0+m1x+m2y+m3x²+m4xy+m5y²とする。
元の全七源は

    {q1,q2}=R+βyq1−αxq2,
    {q1,x}=αx², {q1,y}=1+βxy,
    {q2,x}=1+αxy, {q2,y}=βy², {x,y}=0.

これは公開 `CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json` case7 の
full_six_magnetic_sourceと全係数で一致する。
αβI≠0の稠密主開で

    q1'=(q1+m2/β+m5y/(2β))/I,
    q2'=(q2−m1/α−m3x/(2α)+m4y)/I

とし、tensorをI倍すると、全六括弧がR=1の同じ源へ正確に移る。
全逆写像、Jacobian、全六式を記号m0,…,m5,λのまま検証した。
次に

    U1=q1', U2=x, V1=y, V2=y−q2'

と置く。この実linear変更により定数部分は
`J=∂U1∧∂V1+∂U2∧∂V2` となり、目標の全λ源は

    {U1,U2}=αU2²,
    {U1,V1}=1+βU2V1,
    {U1,V2}=−U2V1−βU1V1−αU2V2,
    {U2,V1}=0,
    {U2,V2}=1+αU2V1,
    {V1,V2}=βV1².                                      (1)

となる。以降の曲線は全λで多項式なので、まず(1)の全λを認識する。
その後に全七源の既約性と閉性を使ってαβI=0も戻す。
磁気の特殊値だけを一般源と呼んでいるのではない。

## 2. 真正な全三parameter親族と実SNC極限

親座標を(u1,u2,v1,v2)、自由係数をA,B,cとする。
次の全三parameter二次Poisson源を取る。

    {u1,u2}=−c u1u2,
    {u1,v1}=1+A u1v1+B u2v2,
    {u1,v2}=c u1v2,
    {u2,v1}=(B−c)u2v1,
    {u2,v2}=1+B u2v2,
    {v1,v2}=(B−c)v1v2.                                 (2)

全四strict Jacobiが自由係数A,B,cのまま恒等零である。
Pfaffianは正確に

    −(1+B u2v2)(1+B u2v2+A u1v1).                       (3)

AB≠0なら斉次化した二因子はrank三とrank五で、相異なる既約quadricである。
ただし最終帰属には全H+2Q分類を使う必要すらない。
座標重み `(h,u1,u2,v1,v2)=(3,0,1,0,1)` による実PGL作用を取ると、
全定数項と交換項B u2v2∂u1∧∂v1だけが正重みで消える。
極限はlogcanonicalで、実五次式は `−AB*h*u1*u2*v1*v2`。
従って五独立超平面の実SNC log-symplectic点である。
connected PGL作用は各全空間成分を保存するので、(2)の点を通る同じ成分がSNC点を含む。
MPS Theorem1.4の全period・全SNC scopeから(2)はMPS40へ入る。
特定の成分の一意性、cの非零性、generic eigenvalue ratioは不要である。

これは `TWO_QUADRICS_SUPPORT.md` §3のB≠0全源の実SNC証明と同じである。
同noteの座標(x,y,z,w)=(u1,v1,u2,v2)で、

    C=D=1, E=0, F=B,
    p=−c, q=c, r=c−B, s=B−c,
    v1=p+q=0, v2=p+r=−B

と直接一致する。親period又はその積を独立Taylor係数として自由化していない。

## 3. 全λの実symplectic confluence

εを曲線parameterとし、(2)の真正係数を

    A=−ε(2ε+1), B=−ε(ε+1), c=εα                       (4)

とする。ε≠0で実symplectic linear変換

    U1=u1+u2/ε, U2=u2,
    V1=v1, V2=v2−v1/ε                                (5)

を作用させる。U側はG=[[1,1/ε],[0,1]]、V側はG⁻ᵗであり、Jを正確に保つ。
全負冪が相殺し、次の文字通りの多項式族になる。

    {U1,U2}=αU2²−εαU1U2,
    {U1,V1}=1−ε(2ε+1)U1V1−ε(ε+1)U2V2+βU2V1,
    {U1,V2}=εαU1V2+(2ε−β)U1V1−αU2V2−U2V1,
    {U2,V1}=ε(β−ε)U2V1,
    {U2,V2}=1−ε(ε+1)U2V2+αU2V1,
    {V1,V2}=(β−ε)V1(V1+εV2).                          (6)

ε=0は全λについて(1)そのものである。
全四strict Jacobi、hによる全二次斉次化、全350 projective Jacobiを、ε,λの記号のまま直接検証した。
実Pfaffianは

    −[(ε+1)U2(V1+εV2)−1]
     *[ε(2ε+1)U1V1+ε(ε+1)U2V2−εU2V1−1].              (7)

ε≠0,−1,−1/2では(4)のABが非零である。
よって(5)の前の親は§2の実SNC極限を持ち、(6)も既知MPS40の閉和に入る。
有限閉和の閉性からε=0も入り、これが全λに成立する。
除外した三つのεは曲線のgeneric親chartだけの除外であり、端点を捨てていない。

全七源はλ,m0,…,m5を持つ既約affine源である。
§1のαβI≠0は非空稠密なので、その全像がMPS閉和へ入ると、
αβI=0を含む元の全源像も閉性で入る。これで全quadratic_AのE義務を閉じる。

## 4. exact証拠と発見経路の区別

`CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_SOURCE_20261107.json` は、元の全七源、全六磁気正規化、
symplectic標準化、五係数semi-simple sourceの全Jacobi idealの両方向certificateを保存する。
その三componentはb=e=0、e=0/f=b、b=0/a=eである。
この全五係数分類や、最初の固定一次eigen-arcの負の結果は、最弱の最終証明には不要である。
上の明示親(2)、実係数(4)、実変換(5)、全族(6)、全源正規化だけで十分である。

`CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_SPECTRAL_ARC_20261107.json` の第二curveが(6)である。
全4/all350、pole相殺、全λ端点を0.871秒で検証した。
最初のarcはeigenvalue中心の二次係数を欠いていた。
実symplectic表示のpole二から、その係数だけが端点を変えると分かり、
一回の有限四unknown系で(r0,r1,c0,c1)=(1,2,0,λ−1/2)を求めた。
これはhigher Jacobi jet、局所Mora、unit消去の再試行ではない。
旧失敗report、旧十八変数idealとその未解決状態は不変である。

独立レビューでは最弱の§§1–3を別算術で再構成すればよい。
形式支持の零次元性又はfullsymplectic source分類を追加条件にしない。
