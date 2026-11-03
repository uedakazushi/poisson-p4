# W1 の actual S61 交叉認識：多項式移送と全完成環の独立検証

2026-11-03。型の一致から交叉を推測せず、元のW=1テンソルを旧universal actual imageへ
実際に移すことを検証する。結論は、全非零Wにおけるactual Dのcone不変量
`(dim,edim,depth,e)=(28,37,21,2)`、主平方零nilradical、21次元の滑らかなembedded支持である。
W=0の元のB18点、環の直和分裂、globalな成分数には拡張しない。

## 物理的なframeとimmutable座標adapter

独立コードはproducerを呼ばず、有理数添字演算で元のW=1点の全150係数を作る。
行列

    (old h,old x,old y,old z,old w)
       ↦(old x,−old h,old y,old z−old x,old w−old h)

の逆行列と行列式1を別に計算し、点(0:1:0:1:0)をe0へ送ることを確認する。
二次単項式の逆座標代入とexterior作用を直接掛けて全150係数を比較する。
得られるテンソルは正確に

    x0²J−x0*x1 ∂2∧∂3+(2x0*x3+x1²+x3²) ∂2∧∂4.

次に25個のEuler-wedge行から垂直24次元空間を新たにRREFし、immutable v02のpivotと一致することを確認する。
その商でperspective pivotは全て零であり、S61の61自由座標を直接抽出できる。
元のencoded eliminated-raw式による150再構成、全134関係を独立に有理数評価する。
S61.normalizeやproducerの再実行に依存した確認ではない。

β16→intrinsic24 adapterを別に掛け、旧sourceのp=s=0のcritical β24に文字どおり一致する。
旧全61source式へ

    p=s=r=0, k0=1, k5=1/3, 他のk=0

を代入してCA/CB/Nの三族で同じ61点になることも確認する。
これはH=x0(x1x4+x3²)+x1²x3+x3³/3の点であり、別のSp4 orbit探索を介さない。

## 新しい点で全actual chartが使えること

旧fiber座標Kを用いたsource→targetの選択行Jacobianを、保存された全61多項式から新たに微分・評価する。
CA/CB/Nの三つ全ての選択行列式は−72である。
βの七群方向の単位minorも別の有理数逆行列で確認する。
このminorはβだけに依存し、同じcritical βでKの値を動かしても失われない。

旧actual image表示では、fiberの全source変数はK,p,sの多項式で戻せる。
全109/109式はregular-sequence Koszulとproductによる**全image kernel**であり、
指定旧点の完成環だけに適合させた近似式ではない。
新点のactual欠損odd行列のrankは7/7/6で、k0,k5方向の列は零である。
そのためこの自由変数の移動後も、a,sを加えた九/八生成regular sequenceの同じ全像表示が使える。
proper有限source chartと七群方向の形式積を、このunit条件によって新点で再適用する。

共通49 targetにおける点は正確に `v4=v7=1`、他は零となる。
保存された全218式を別の有限多項式演算で評価し、全218×49 Jacobianの全係数を再構成する。
CA、CB、和idealのrankは28、28、29、nilpotent Fの微分を加えたrankは30である。
v4,v7をそれぞれ1だけ移した全218多項式と、保存されたtranslated有限表示も係数ごとに一致する。

## globalな多項式等式の局所化

既受理の独立M2 gateは、QQ[49変数]そのものの上で

    F²∈ID, ID+(F)=IN, ID:F=Qemb

を検査している。今回のコードはそのscriptが元の多項式環で定義されていることと全入力hashを確認する。
新しいM2を実行する必要はない。元の点の形式的な冪級数恒等式を、離れた点へ勝手に平行移動していない。

Qembはv1,v4,v5,v7以外の45独立一次式のprimeである。
新点で全て零となり、45次元の独立性も別に確認する。
したがって局所化後も `(F)≅(R/Qemb)F` は非零である。
F²=0とprimeな全N像 `ID+(F)=IN` から、nilradicalは正確に主平方零ideal(F)となる。
Qembの自由v4,v7方向の移動はこのprimeを変えない。

## 全環、depth、associated primesとsmooth因子

critical sliceでNのregular odd idealは八生成である。
元の全kernel議論を新しいsource点へ適用すると、NはAc7に三つのsmooth自由方向を加えた型となり、
dim11、depth4、多重度2を持つ。nilradical加群は正則4次元の `R/Qemb` と等しい。
従って

    0→R/Qemb→R/ID→R/IN→0

でdepth lemmaはdepth≥4を与える。
local cohomologyではH³(R/IN)=0、H⁴(R/Qemb)≠0がH⁴(R/ID)へ単射するため、depthは正確に4である。
associated primesはNの極小primeとQembの正確に二つとなる。
前者は唯一の極小prime、後者は実際の部分加群のannihilatorであり、
完全列によるAssの包含から他のprimeは生じない。
QembはNより低次元なので、top multiplicityはNの2を保つ。

七つのβ-orbit方向に加え、symplectic値の六方向とperspective四方向を戻す。
valueはe0で文字どおりJであり、この受理済みの十方向のsmooth分離のunit条件を満たす。
coneでは合計17変数、projectiveではscaleを除いて16変数を加える。

* cone Dはdim28、edim37、depth21、e2。
* projective Dはdim27、edim36、depth20、e2。
* nilradicalは非零で平方零、その全加群は正則21次元環上free rank1。
* Assは極小Nと滑らかな21次元embedded prime（projectiveでは20次元）の二つ。

全完成環は、固定した49変数218式を上記の新点へ移した商の完成に17smooth変数を加えたものとして
実際に指定される。trivial-extension直和への分裂は主張しない。
schemeを保つ全非零WのGL/scalar群恒等式で、この全結果をW≠0へ移せる。
元のrank2 B18点W=0はS61のこのrank4 chartの外であり、今回のdepth21を移せない。

専用コード：`code/ambient_W1_S61_intersection_recognition_review_20261103_v63.py`。
flat receipt：`reports/AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63_REVIEW.json`。
60秒・2GB上限、source jet・M2・GBの再計算なし。
