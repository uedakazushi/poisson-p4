# nil22 quadratic_A：同一既存軌道の完全18変数横断模型

2026年11月7日、`E_multiple_hyperplane`。全七変数源の成分帰属は未解決。
本ノートは切片の厳密な有限変数削減を与えるが、被約支持を判定したものではない。

## 1. 既存anchorを再利用する

全源の七変数は六磁気とλである。`NIL22_IMAGE_SOURCES_AND_COMPONENTS.md` §10 の実fiber平行移動は
全六磁気を不変量 `I=m0+m4` へ帰着する。既存証明書のanchorは
`λ=7,m=(1,2,3,5,7,11)` なので I=8≠0。
共通fiber/tensor scaleにより、**同じPGL/scalar軌道**を `λ=7,F0=1` と表示できる。
新たな有理点を探索していない。

この表示では `(h,q1,q2,x,y)` の重み `(0,1,−1,1,−1)` に対して tensor重み0である。
全150係数と全350 projective Jacobiを使い、Euler24、実GL5、scalarと**全七source微分**を保持すると

    rankD=101, rankG=47, rankF=48

となる。真正法線は一つ、そのトーラス重みは0である。
従って正重みだけの有限IFT次数評価は適用できない。

## 2. 完成切片をトーラス固定部分へ厳密に縮める

Dの独立101標準座標列と、weight0の真正法線Nで102次元の線形空間Lを取る。
`rank[F,L]=150` を確認したので、Lは48本の実線形cutに対する真正切片である。
各座標列はトーラスの固有ベクトルであり、Lはトーラス不変となる。

独立101 Jacobi成分を implicit 方程式として選ぶ。各成分も固有weightを持つ。
implicit Jacobianの101×101小行列は可逆なので、完成局所環では唯一の形式写像

    u=u(n), u(0)=0

で implicit変数を解ける。normal n のweightは0だから、n に対するトーラス作用は自明である。
implicit方程式の同変性と解の一意性により `t·u(n)=u(n)` が全tで成立する。
従って**非零weightを持つ84個のimplicit変数の形式解は恒等的に零**である。
これはnormalの高次項の切断ではなく、完成環における厳密な消去である。

残るweight0 ambient係数は34個で、同blockのD階数は17、F階数は16。
真正固定切片は `17 implicit +1 normal =18` 変数となる。
この18変数の二次斉次二ベクトルを元の全350 Jacobiへ代入したところ、
非零の方程式は45本であった。他の305式は恒等的に零。
従って完成した**全102切片**の支持判定は、完全な18変数45式の局所idealへ帰着する。

公開入力は
`reports/CLASSIFICATION_NIL22_QUADRATIC_A_WEIGHT_FRAME_20261107.json` と
`reports/CLASSIFICATION_NIL22_QUADRATIC_A_FIXED_MODEL_20261107.json`。
全350の制限式、45生成元、17個のGauss正規化式、全minor逆行列、native入力を保存した。
17個のGauss式の一次斉次部分は正確に各implicit変数 z_i である。
もし元の45生成元idealの一員で最小斉次部分が `c n^k,c≠0` のものを一つ証明できれば、
全18初期式の根基は原点となり、固定切片・全切片とも局所次元0となる。
48本のcutのKrull不等式と実全source F48が一致し、全七源の閉包が全周囲23成分であると結論できる。
この最後の一員は現時点では未取得である。

## 3. 実行上の結果と停止した経路

有限契約 `jobs/nil22_quadratic_a_fixed_support.json` のnative localMora `liftstd` は
既存Singularを用い、18変数45式に対し85秒で未完了となった。
支持次元、正規冪、成分性について**何の結論も出ていない**。
実行inputとproducer bytes/hash、全数学模型を上記public JSONへ保持し、
`CLASSIFICATION_NIL22_QUADRATIC_A_FIXED_RESOURCE_LIMIT_20261107.json` に上限を明記した。

最初のloader実行は既設shared-library環境未activationで起動せず、
既存cloud-setupのLD_LIBRARY_PATHを明示して修正した。再セットアップは行っていない。
その後に実際に走った85秒jobが上の未完了結果である。

このlocalMora経路は一回の上限で停止する。より長いjet、別anchorの点列、予算増加は自動追加しない。
変換行列を持たない完全stdや別の厳密構造定理は別算法だが、現在の支持結論とは区別する。
未解決の全λ源に対する実族・既存成分認識は引き続き必要である。

## 4. 厳密なunit消去とliteral-linear模型

nativeMoraの再試行は行わず、実際の生成元にあるunit一次係数による
局所環の有理同型を使う別算法を一回試した。
`classification_local_unit_elimination.py` は、変数zに一次の生成元 `u z+v` で
`u(0)≠0` のものだけを使い、全45式を保って置換する。jetではない。
しかし一般的な疎性基準の逐次消去は有理式が膨張し、90秒で未完了となった。
完成した消去鎖・最終一変数ideal・支持証明は得ていない。
実行bytes/hashと上限は `CLASSIFICATION_NIL22_QUADRATIC_A_UNIT_RESOURCE_LIMIT_20261107.json` に保存した。

45式中のliteral一次式は12本あるが、**そのrankは7**である。
一括の定数Gaussian消去は厳密に完了し、18変数から11変数の全45式模型へ移した。
`CLASSIFICATION_NIL22_QUADRATIC_A_LITERAL_LINEAR_MODEL_20261107.json` は全係数行列、
七変数の逆置換、残る全45式を保存する。これは全idealの実同型でありパラメータ開を落とさない。
元の疎な式には `z2=n z1−43n²`、`z7=−2n z11−56n²` があるが、
それらだけから全implicit変数が有理的に消せるとは結論していない。

二つの局所支持算法が有限上限で閉じなかったため、この経路は停止する。
消去順序の列、別anchor、より長いnormal jetを自動追加しない。
成分性と全λ源の帰属は未解決として保持する。
