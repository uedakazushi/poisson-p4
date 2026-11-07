# nil22 image double_B の全六磁気源は CPLUS_EXT に含まれる

2026年11月7日。`E_multiple_hyperplane`。全 Jacobi/Pfaffian の直接計算と別作者の
Fraction-only レビューが PASS した。新成分は追加しない。

`NIL22_IMAGE_SOURCES_AND_COMPONENTS.md` §3 の double_B を扱う。
座標は `(q1,q2,x,y)` とし、任意の六係数

    R=m0+m1x+m2y+m3x²+m4xy+m5y²

を保持する。元の源の二anchorは

    U0=(−x²/4−1, −xy/2+2x),
    V0=(−xy/4+x, 8−y²/2),
    {q1,q2}=R+(1−y/4)q1.

以下の実族は全源を一度に送る。ε=2を除き

    a=(ε+2)/(ε−2),
    Uε=(−x²/4+2/(ε−2), −xy/2−2ax),
    Vε=(−xy/4−ax, 8−y²/2+8εx²/(ε−2)),
    Aε=y−4+εx²,       B=y+4+2x²

と置く。

    {q1,(x,y)}=Uε,   {q2,(x,y)}=Vε,
    {q1,q2}=R+(−y/4−a)q1,   {x,y}=0.                    (1)

直接に

    [Uε,Vε]=(−y/4−a)Uε,
    det(Uε,Vε)=AεB/(2−ε)

を満たすので、(1) は全 ε と全 m に対して Poisson である。
別の全350 projective Jacobi も確認した。
ε=0では Uε=U0,Vε=V0 であり、全六磁気の元の源へ**全六括弧が正確に**戻る。
唯一の族の分母 ε−2 は0で非零であり、磁気係数の反転又は特殊化を一切行わない。

h により斉次化した実反標準因子は

    Dε=−h QA QB/(2−ε),
    QA=hy−4h²+εx²,   QB=hy+4h²+2x².                    (2)

ε≠0,2 では QA,QB は相異なるrank3二次式で、h=0を同一点 `[0:0:1]` で共通接線に持つ。
変換 `(h,x,y)=(X,Y,−2Z−4X)` によって

    QB=−2(XZ−Y²),
    QA=−ε{(2/ε)XZ−Y²+(8/ε)X²}.

従って `shared_tangent_generic` の全配置に入り、配置変数は
`u=2/ε≠0,1`、非零X²係数は `8/ε` である。
このX²係数の規格化に平方根を使っても複素数上の可逆変換で足り、族の正則性や端点には影響しない。
`TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md` §1,5 の全free-frame・全六磁気源定理により、
ε≠0,2 の全 (1) は既存24次元成分 `R4_CPLUS_EXT` に属する。
これは配置だけを合わせた一点の帰属ではなく、全Poisson源に適用する既証明である。
u=1/2に対応するε=4も、その既受理の別frameで含まれる。

各固定mに対し CPLUS_EXT の閉性を ε=0 へ適用すれば、元のdouble_B点も同成分に入る。
従って**全六磁気パラメータと全境界を含む double_B 源の全支持**が被覆された。
磁気正規形の二つの不変量を選んだ部分だけの主張ではない。

公開のproducerは `code/classification_nil22_double_b_cplus_limit.py`、
結果は `reports/CLASSIFICATION_NIL22_DOUBLE_B_CPLUS_LIMIT_20261107.json`、
有限契約は `jobs/nil22_double_b_cplus_limit.json`（直接計算0.232秒）。
独立レビューは `code/classification_nil22_double_b_cplus_limit_review.py`、
`NIL22_DOUBLE_B_CPLUS_LIMIT_REVIEW.md` と対応report/jobに保存される。
別のFraction算術で全4/350、全六源端点、実Pf、両rank3行列、共通接線配置、
全Saito-frameと磁気の可逆6×6変換（行列式16）を再構成した。

磁気ゼロ点の graded frame では真正法線に重み0が残るため、有限positive IFTは使っていない。
被覆は式(1)という実際の有理Poisson族と、既存の全源包含・閉性のみから従う。
