# nil22 double型：全六磁気の実affine正規形

2026-11-07。`E_multiple_hyperplane`。
残る E 認識のため、double_A/B/C の全六磁気係数に対する fiber affine 作用を有限に解いた。
独立 code は `code/classification_nil22_image_recognition.py`、報告は
`CLASSIFICATION_NIL22_IMAGE_MAGNETIC_NORMALFORMS_20261107.json`、
契約は `jobs/nil22_image_magnetic_normalforms.json`。
全係数の exact 演算が0.302秒で PASS。既済有限分類・成分証明書は変更していない。

P を `F∂q1∧∂q2+ΣMij∂qi∧∂pj` とし、二行 U,V、F の q 係数を α,β とする。
任意の二つの affine 式 L1(p),L2(p) で q'=q+L(p) とすると、M はそのままで磁気項だけが

    F0'=F0−αL1−βL2−V(L1)+U(L2)

へ移る。これは determinant1 の実 affine PGL 変換である。
六未知係数に関する全6×6行列を作り、核・余核・全元係数に対する逆変換を解いた。
その逆変換を実際に tensor へ代入し、全六括弧の一致も検証した。

作用階数は double_A/B/C の順に5,4,3。
各残磁気不変量と正確な L1,L2 は報告書に保持する。
これらは source の一部分を選んだ結果ではなく、全六係数の恒等式である。

double_A の唯一の不変量は

    I=m1−m3−4m4/3.

全六源は追加の係数反転なしで F0=I x へ移る。
I非零なら q'=q/I と tensor の I 倍により F0=x へ正確に移る。
I=0 の全補集合は F0=0 へ移り、F0=τx の全 Poisson 多項式曲線の τ=0 endpoint となる。
実 Pf が磁気に依存しないので、この曲線で零 tensor の問題もない。

従って **double_A 全六源の projective PGL 飽和閉包は、F0=x の一つの実軌道の閉包に等しい。**
全点が同一 actual orbit に入るとは主張しない。
既存 witness の orbit23 ともこの正規化で対応するため、全帰属はこの canonical orbit の認識へ還元できる。
現時点では ambient componenthood や既知成分への membership は未解決である。

rank6なら全磁気を0へ戻して generic q-scaling torus を使う予定だったが、実階数は5。
q の尺度は I を動かすので、F0=x の generic point にその torus はない。
F0=0 境界だけの torus から generic 点の有限 weighted model を断定しない。
二次障害零を理由とする高次 jet 反復には進んでいない。
root の優先指示により、この認識作業はいったん保存し、nil31-long 第一源の成分監査へ移る。

## 全affine N∧Y pencil の有限診断

第一 nil31-long 成分監査完了後、canonical double_A の F0=x に戻り、
全20係数の affine Y に対して `P+τ(N∧Y)` の一次 compatibility を完全に解いた。
`CLASSIFICATION_NIL22_IMAGE_DOUBLE_A_TWISTS_20261107.json` と
`jobs/nil22_image_double_a_twists.json` に全核・self-Poisson条件・実pencilを保存した。
1.791秒で PASS。

compatibility kernel は正確に六次元で、

    Y=L1(x,y)∂q1+L2(x,y)∂q2,

L1,L2任意affine式だけである。全self-Poisson条件も恒等零。
この場合 N∧Y は `(xL2−yL1)∂q1∧∂q2` だけなので、全pencilは単に磁気を動かす。
全周囲微分に対する真の法線像も零で、cofactor は全τで同じ N のまま。
従ってこの **全affine Y** の経路は既に正規化した同じsource内から出ず、
半単純場への実曲線を与えない。SS20曲線数は0で、成分帰属の結論は増えない。

これは他の全Poisson族の不存在ではない。
例えば多項式 W が `N∧W=0` を満たしても、`W=N∧Y` に多項式Yが存在するとは限らない。
N の成分が生成する ideal は(x,y)で、q方向の二次式をqq括弧へ加える cohomology が残る。
この欠落は無限jetではなく、有限Koszul表示の追加六クラスに対応する。

## 欠落六クラスを加えた全二次N∧W=0も磁気だけ

この欠落を明示的に補い、全二次 W に `N∧W=0` を課した全空間を

    Wqq=任意の四変数二次以下（15係数）,
    Wqi,pj=pi Rj(q,p), Rj=任意affine（各5係数）,
    Wpp=0

という全25係数で保持した。x,y が互いに素なので、N∧W=0 の二つの cross 関係から
上の表示が必要となり、逆も恒等的に成立する。従って追加 slice 仮定ではない。
`jobs/nil22_image_double_a_koszul_pencils.json` の有限契約で全 compatibility を解いた。
`CLASSIFICATION_NIL22_IMAGE_DOUBLE_A_KOSZUL_PENCILS_20261107.json` は2.464秒で PASS。

核は正確に六次元で、二次式の係数 indices は `{0,3,4,12,13,14}`。
すなわち

    W=R(x,y)∂q1∧∂q2,  deg R≤2

だけである。全 cross 係数と q 依存の cohomology は零に強制される。
この全核の self-Poisson 条件も恒等零で、全pencilは元の全六磁気 source 内にとどまる。
従って **N∧W=0 を満たす全二次 Poisson pencil** の経路からは、別の cofactor や既済SS20源へ出られない。
法線像零という接空間情報だけでなく、完全な25係数表示とその全核がこの結論を与える。

これは二回目の非閉鎖経路の有限診断なので、この N-wedge route をここで終了する。
N∧W≠0 の別の実Poisson族や他の幾何認識は排除していない。
double_A の既知成分帰属は引き続き未解決で、追加jet・無制限pencil係数scanは行わない。

## T39/T40全graphの全source検査

別経路として、Aは canonical F0=x、Bは残る全二変数 `F0=I0+I1 x`、
Cは残る全三変数 `F0=I0+I1 x+I2 x²` を保持して faithful S61 へ移した。
原点のPfは各sourceの全磁気に依存しない非零定数。
固定したcrossing成分をpivotとするDarboux frameとperspectiveは、全残磁気に関して
定数以外の分母を持たない。このため係数零境界を含めて同じchartを使える。
全150再構成、全61座標、全134式、全350Jacobi係数を確認した上で、
T39/T40の全45graph式を厳密に評価した。

`CLASSIFICATION_NIL22_IMAGE_FULL_GRAPH_RECOGNITION_20261107.json`、
`jobs/nil22_image_full_graph_recognition.json` に全frame・perspective・61座標と全残差を保存した。
5.744秒で PASS。最短の定数残差は次の通り。

| canonical source | T39残差 | T40残差 |
|---|---|---|
| A:F0=x | u2で −9/1024 | u2で 5/256 |
| B:全I0,I1 | u2で 3/16 | u2で 1/4 |
| C:全I0,I1,I2 | u6で −1 | u6で −1 |

従ってB/C全sourceは両成分に入らない。AではI非零の全actual orbitが両成分に入らない。
AのI=0境界について、この負結果だけから非所属とはしない。
これらは二成分に対する完全なgraph検査であり、他MPSへの非所属や新成分性を意味しない。
N-wedge route は終了したままで、認識を全projective族・P5境界の幾何へ切り替える。

## double_A の既知97内の候補はP5三閉包だけ

`CLASSIFICATION_NIL22_IMAGE_DOUBLE_A_KNOWN_FRONTIER_20261107.json` と
`jobs/nil22_image_double_a_known_frontier.json` は0.581秒で PASS。
ここでは各既知成分への帰属候補を有限に絞るだけで、新成分やmembershipを主張しない。

canonical double_A の実軌道は23次元であることを、全GL5・Euler・scale行列の階数48から再現した。
その三次因子 C の y 二次式としての判別式は

    (4x−1)(4x+3)³/256.

C(x,y) は絶対既約で、唯一のaffine特異点は `(x,y)=(−3/4,−4)`、Hessian階数1のordinary cusp。
実際、X=x+3/4、Z=y+4+(32/3)X とすると
`C=−3Z²/128+X²Z−(32/3)X³`。平方完成後の三次接項が非零なのでnodeではない。
既済nil22新23二つはnodeを持つcubicと23次元稠密軌道を持つ。
同じ23成分なら二つの稠密開軌道が交わるため、cusp/nodeの相違によりAとは異なる。

次にMPSの各24次元成分について、全period/smoothing sourceを保持したまま
対角PGL作用の全係数恒等式を検査した。
各smoothable edge i,j の重みは `ei+ej−Σθk ek`。
これらを全て消す対角場は、全smoothing completionの各monomialを厳密に固定する。
global scalarを除いた対角stabilizerの次元下界は4−rank(edge weights)。
完全なsourceでの恒等式から、次の全成分上の軌道次元上界が従う。

| MPS成分 | edge weight階数 | 対角stabilizer下界 | 全成分上のPGL軌道上界 |
|---|---:|---:|---:|
| T02,T03 | 1 | 3 | 21 |
| T16,T17,T21 | 2 | 2 | 22 |
| T36,T37,T38 | 4 | 0 | 24 |

T17/T36/T37/T38では既存の厳密な高次smoothing completionを用い、未完成の最初の和で代用していない。
他の四源は既済の全generic sourceを用いた。
sourceの稠密性・成分ラベルは既存MPS証明を依存とする。
群作用の軌道階数上界は閉条件なので、generic source上の上界はその全成分閉包へ延びる。
従って23次元軌道を持つAは、最初の五成分に入らない。

既存96と追加済みnil31-long22を合わせ、全97 IDごとに次を監査した。

* rank2は閉rank条件、次元≤22は23次元軌道を含めないことにより排除。
* 23次元の既約成分に入るならAの23次元軌道閉包と等しい。
  Aの全閉包は閉条件χ=H²C3に入るが、既存23の大半は被約/SNC点を持つため等しくない。
  zero-logは閉条件χ=H⁵で排除。node cubicの二つは前述のdense orbitの相違で排除。
* T01は全閉包の五一次因子条件で、T39/T40は全graphの非零残差で排除。
* 次元≥24の非MPS全成分は、その監査済み一般因子型から、全閉包が次数≤2の因子の
  proper projective積写像の像に入る。Aの既約三次因子はそこに入らない。
* 残るMPS24五成分は上の全軌道階数上界で排除。

従って **既知97のうち canonical double_A を含み得るのは T36,T37,T38 の三つだけ**。
全六源はcanonical軌道閉包なので、その既知閉和への全帰属もこの一点の三閉包認識へ正確に帰着する。
三候補のどれかが実際に含むという証明はまだなく、未知成分の存在や次元も結論していない。
今後必要なのは、P5の既済全projective sourceからの実PGL退化、又は元150係数のh行も動かす
明示的な全Poisson族である。停止済みのaffine N-wedge検査や高次jetの延長には戻らない。
