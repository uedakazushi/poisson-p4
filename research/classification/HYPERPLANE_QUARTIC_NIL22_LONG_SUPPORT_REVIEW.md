# nil22・長translationの全源次元上界：別作者独立確認

2026-11-07。`HYPERPLANE_QUARTIC_NIL22_LONG_SUPPORT.md` の全20係数源を独立再構成した。
結論は固定cofactor源の次元≤9、全射影飽和像の次元≤24であり、
既約quartic条件を課す前の全支持、可約・零Pfaffianの全境界を含む。
Eの全成分帰属は主張しない。

独立実装は `code/classification_nil22_long_support_review.py`。
producer、CAS、Singular、NumPy、SymPyをimportせず、既存の独立Fraction疎多項式parserと
行列演算だけを再利用した。実行前gateは `jobs/hyperplane_quartic_nil22_long_support_review.json`。
実行3.870秒、最大RSS25,204KiBでPASS。
証明書は `CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_SUPPORT_REVIEW_20261107.json`。
reviewerの完全bytesも証明書内に保存した。

全90係数線形源自体には既済の独立110行再構成を入力とした。
その全20parameter括弧から四つの元Jacobiを改めて生成し、全係数と七つのprimitive core生成元を照合した。
従ってcoreの式だけを外部の理想として受け入れたものではない。
core・三criticalcore・rankA≤1閉層の五証明書では、
元生成元、全670個のideal membership lift、全127個のlift先polynomialを再構成した。
各恒等式 `G_j=Σ_i L_ij I_i` はQQ上の疎多項式として一致した。
grevlexのleading monomialも各polynomialそのものから独立に取り直した。

表示されたGが入力理想Iに属するので、全LM(G)が生成するmonomial idealはin(I)に含まれる。
従ってそのmonomial idealの次元がV(I)の次元の上界である。
G全体が完全なGroebner basisであることやSingularのdim出力は必要ない。
八core変数の全有限subsetを小さいcardinalityから列挙してminimum vertex coverを求め、
全core、q=−1、q=−1/2、q=−1/3、q=−1/2かつrankA≤1の次元上界
`2,1,2,1,1` を確認した。反転するcore元は無い。

全元Jacobi係数から2×5行列A、5×7行列Bと全右辺を取り直した。
Aのrank≥1は定数恒等式−4で全点に成立する。
Bの二行四列block Cが零なら六core変数が一意に強制されることを非零定数6次minorで確認し、
その全代入で元Jacobi係数第25番が実際に3となることを検査した。
従って全core点でC≠0である。
regular三pivotと三critical二pivotの各determinantを再計算し、
Cの二行はそれら全pivot列で零であることも確認した。
よって全aと全coreでBのrankはregular≥4、critical≥3。
これは一般点又は標本でのrank下界ではない。

以上を全五層へ適用すると、各層の次元上界はそれぞれ
`2+4+3`, `1+4+4`, `1+4+4`, `2+3+4`, `1+4+4` で、全て9である。
rhsの整合条件がfiberを空にする場合も上界は増えない。
rankA=0は既に不可能、rankAは高々2なので、全点を覆う。

全affine作用の20×20行列も列ごとに独立構成し、線形共役部分rank8、全rank12を確認した。
cofactor Z=y∂x+∂y+w∂zについて
`diag(λ²,λ,λ,1)` がZをλZへ送るため、全飽和源はtensorの非零scaleを保つ。
従って超平面選択四次元とprojective scalarを数えて
`9+12+4−1=24`。
未知W優越成分は次元25–27なので、この全型はその一般支持にならない。
同じ全source上界はH+二既約quadricにもそのまま適用できる。

原noteの二段fiber、全critical層、initial idealへの包含、affine軌道と射影化の各橋も
独立に読解確認した。失敗した全20変数GB試行には依存せず、GBを再実行していない。
