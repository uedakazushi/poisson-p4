# Ga9非零定数cの全十源は新しい23次元成分を与える

2026-11-07。`E_reduced_nonnormal`。
`HYPERPLANE_QUARTIC_NIL211_SHORT_TRANSLATION.md` §7の**全十parameter源**の
PGL飽和被約像閉包は、全Poisson空間の23次元既約成分である。
根拠は、既存b=2点での完全なQQ係数C*切片と四つの明示ideal-membership恒等式。
任意jetからの推測、QQ(b)で未完了の計算、特定のsourceだけの次元上界は使わない。
既存97成分との相違は§5。別作者Fraction reviewはPASSし、正式な計数は統合担当が行う。
Ga9の他四源、全H+quartic又は全Poisson分類の完成は主張しない。

## 1. 全源、rank4条件と既約性

固定modular場は `Z=y∂x+∂z`、u=x−yzとする。
受理済みの全源はr≠0と九自由変数b,α,p,q,A,D,G,H,kからなり、

    B=b−1/4,
    F2=−(1−b)B/(2r),
    E=α−[(1−3b)A+qB]/r,
    φ=r y⁴/4+A y³+B y²w+D y²+E yw+F2 w²+G y+H w,
    Φ=u+φ,
    Vy=−α−by, Vw=p+qy+bw+r y²,
    A2=φ+Vy φy+Vw φw−k, Vu=−u−A2.

元四変数への逆変換後の全六括弧が二次以下、全Jacobiが恒等零であること、
元28係数源の全支持がこの十源に尽くされることは
`NIL211_SHORT_TRANSLATION_REVIEW.md` §4で独立確認済み。
このsourceはr≠0の既約parameter空間であり、任意のb又は低次parameterを反転していない。
実affine PfaffianはF=k−u−φ。
斉次化F4のx係数は−h³、定数項はr y⁴/4を持つためhで割れない。
Gauss補題でF4は全parameter域で絶対既約であり、D=hF4は相異なる二因子H+F4である。

## 2. SAME既存anchorと真の全周囲切片

r=1,b=2、他parameter=0を取る。既存k=1固定点からx'=x−hでk=0へ移しただけであり、
新しい有理点を探したのではない。8倍した全h-Casimir括弧は

    Pxy=14hw−2y²,
    Pxz=8hx+35w²+8yz,
    Pxw=8hz−12wy,
    Pyz=16hy, Pyw=8h², Pzw=16hw+8y².

全150係数、全350 projective Jacobiを保持すると

    rankD=98, rankEuler=24, rankG=47, rankF=48, DF=0.

Fには**全十源の実微分**を加えた。genericなG+∂bPの特殊化だけでFを代用していない。
全kernelは52次元であり、実Fに補な四真正法線N0,…,N3をQQ上で選ぶ。
98本の独立なD標準座標列Uに足したL=[U,N]は102次元で、rank[F,L]=150。
従ってLは実際の48本の独立線形cutの共通零点である。
選択した98Jacobi行とUの98×98 minorは非零。
全block、四normal、全U/row添字、全150directsumはanchor checkerで再構成した。

## 3. 先験的なC*有限性と全350式

weights(h,x,y,z,w)=(0,4,1,3,2)でPのweightは−3。
対角PGL作用とtensor scaleを合わせてPを固定したC*作用を用いる。
相対weightの符号を一斉に反転して、四normal parameterのweightを(1,2,2,3)とする。
全150係数のimplicit相対weightは−4,…,10、全350式は−4,…,16。
F,U,Nと98選択行をweightstableに選べる。

選択98式の形式IFTは一意なのでC* equivariantである。
normalweightが全て正のため、負又は零weightのimplicit座標は原点値零から恒等零。
weight d>0の座標は同じweight dのnormal単項式の有限和となる。
ambientの最大weight10により、IFTは**重み10以下の有限多項式として全ordersで完成する**。
これは計算後に見つけた打切り次数ではなく、実行前のrepresentationからの上限である。

QQ上でこの全多項式を構成し、非零implicit78個を保存した。
全98選択式は多項式として恒等零。
残りもprior最大weight16までの全350行を抽出し、非零179式を得た。
したがってこれは真の102次元切片の**完全な形式ideal**であり、普通の有限jetではない。
日本語証明の有限性補題は `GA9_CONSTANT_WEIGHTED_SOURCE.md` §3と同じだが、
算術を最弱の既存anchor QQへ変えたためQQ(b)での巨大なgcdやΔ累乗は不要になった。

## 4. 四つの多項式恒等式だけで原点支持を証明する

上の全179式をf0,…,f178とする。
三つのpure powerと一つの多項式G0について、QQ[t0,t1,t2,t3]内の明示的な係数を保存した。

    G0 = Σ_j L0j f_j,
    t1⁴= Σ_j L1j f_j,
    t2⁴= Σ_j L2j f_j,
    t3³= Σ_j L3j f_j.

四liftの非零係数数はそれぞれ2,11,11,11。
使う元式は全179式中のindex
`15,16,36,37,38,39,62,63,64,65,66` の十一個だけである。
全lift恒等式を検査し、さらに各係数を適切なweighted homogeneous成分へ射影した後にも再検査した。
G0のleading monomialはt0⁵、他三式のleading monomialはt1⁴,t2⁴,t3³なので初期idealは零次元。
より直接には

    G0(t0,0,0,0)=C t0⁵,
    C=514114689789605014659184048028741872800 ≠0.

三つのpure powerが根基でt1=t2=t3=0を強制し、G0がt0=0を強制する。
原点が全式の零点であることと合わせ、真の切片の根基は(t0,t1,t2,t3)である。
精密なnil指数や完全な非被約環を求める必要はない。

48本の実線形cutのKrull次元不等式から、全150係数の局所次元は48以下。
全実F48による下界と一致する。
全十源×PGL×Euler×scalarの像閉包は既約で、このanchorを通り次元48を持つ。
それを含む全成分も48次元以下なので、その像閉包自身が全空間成分となる。
Euler24とprojective scalar1を引いた次元は23。
この論証は全b、元parameter域全体とその有限閉包境界を含み、generic QQ(b)支持を必要としない。
射影接次元27と成分次元23を区別し、局所一意性・滑らかさは主張しない。

## 5. 既存97成分との相違

全prior97の次元・rank・generic divisor・SNC証拠を
`CLASSIFICATION_GA9_CONSTANT_PRIOR97_20261107.json` に固定した。
成分性を得た後、次元23以外とは次元から異なる。一般rank2成分とも閉rank条件で異なる。

全源でDは平面L={h=y=0}のidealの三乗へ入る。
PGL飽和後の条件「あるP2に沿ってDが三次以上消える」はGrassmannian上の閉incidenceの像である。
従って全成分に成立する。SNC因子はcodimension2一般点で高々二枝しか交わらず、この条件を満たさない。
よってSNC点を持つ全MPS40成分とは異なる。

残る同次元成分のうち、5H、三H+Q、Q+C3、二H+C3、二直線+C3、H+二conicなどは、
対応する掛け算写像の閉像に全generic divisorが入り、全成分のrank4部分にもそのfactorizationが保たれる。
ここには相異なるHと絶対既約F4からなる本源の点が入らない。

同じ1+4因子型を持つ旧四成分はMINIMAL_M、A0、DPM、H_QUARTIC_ESSENTIAL5。
本源の唯一Hのcanonical cofactorは、非零scalarを除いて

    Nhom=y∂x+h∂z,
    Jordan型(2,2,1),  im(Nhom*)=span(h,y),  h∈im(Nhom*).

旧代表はいずれもP=h²J+Q（Qはpure homogeneous quadratic）、h-Casimirであり、
そのaffine divergenceは定数項のない線形場である。
従ってhomogeneous cofactorのdual imageはspan(x,y,z,w)へ入り、h∉im(Nold*)。
MINIMAL_Mはさらにrank1なのでJordan型でも異なる。
独立担当は旧原式から直接発散を取り、ESSENTIAL5/A0/DPMでは全て
`z∂y+x∂w`、MINIMAL_Mでは非零scalarを除いて`x∂y`となることを確認した。
DPMの平方根parameterは形式的に保持し、全て相殺することも検査した。

このflag比較はgenericな全source上で行う。
同じ23次元成分なら両constructible稠密PGL像が共通の非空開を持つが、
uniqueHのdual-image flagが逆なので交わらない。
flag条件そのものを閉条件と仮定してboundary全点へ延長した議論ではない。
以上で既存97とは相異なる新23次元成分を得る。

## 6. 証明書とreview状態

- 全QQ anchor frame：`code/classification_ga9_constant_anchor_frame.py`
- 完全有限IFT：`code/classification_ga9_constant_anchor_support.py`
- 全78implicit・全179残式：`CLASSIFICATION_GA9_CONSTANT_ANCHOR_SUPPORT_20261107.json`
- 最弱四lift：`code/classification_ga9_constant_anchor_certificate.py`, `CLASSIFICATION_GA9_CONSTANT_ANCHOR_CERTIFICATE_20261107.json`
- 契約：`jobs/ga9_constant_anchor_support.json`, `jobs/ga9_constant_anchor_certificate.json`

完全計算1.877秒、四lift0.135秒。既設Singularと既存環境を再利用した。
独立作者は全十源をpotentialから再構成し、全350のexact polynomial、全F/truecut、四lift、
全97との相違を標準ライブラリFractionだけで確認した（4.304秒、RSS23372KiB）。
`GA9_CONSTANT_ANCHOR_REVIEW.md`、`CLASSIFICATION_GA9_CONSTANT_ANCHOR_REVIEW_20261107.json`、
`code/classification_ga9_constant_anchor_review.py` を参照する。
主原稿・global台帳・gitの統合と新成分計数は統合担当が行う。
