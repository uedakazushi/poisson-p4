# 第52 rank4成分の全61座標 dense map：独立監査

2026-11-03。`Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63` を承認する。
18入力・6079nodeの有限有理 circuit は新成分の全S61切断にdenseな像を持ち、
全u1,...,u61を既存規約で出力する。一次像次元15、像上の一般fiber次元1である。
全prime ideal、一次像hypersurfaceの展開式、一般一次部分の五点尽くしは主張しない。

全 canonical150族から一般点へのtranslationを別に組み立て、既に独立監査された汎用normalizerへ渡した。
生成された全6079node、61output index、一般normalization matrix・value pivot・Pfaffian・Sp4 metadataが
保存circuitと完全に一致した。したがってgenericなPoisson恒等式は、実全sourceの全次数恒等式と
正確な射影変換・Euler quotient・S61 gauge定理による。標本の零性から一般恒等式を外挿していない。
全Sp4 big cell十変数、一般評価点四変数、tensor scalarと三族変数が保存されている。

各nodeを独立のdense18次元Fraction dual arithmeticで評価した。
積の微分には逐次積の規則を使い、producerのprefix/suffix実装を再使用していない。
全61値と全61×18微分が保存値に一致し、16次・15次minorを厳密QQで確認した。
全inverseの分母がこの点で非零であることも直接検査し、value pivotは1/36、Pfaffianは−16/81である。

実装のもう一つの照合として、元の物理族のn111 tensorに、source点translationと評価された
normalization matrixの積を直接pushforwardした。独立のmodel.normalizeで得る全61値は一致し、
元の134式、全150再構成係数、直接の射影Jacobiも全て成功した。
genericな像は既知25次元成分の16次元S61切断に含まれ、微分rank16の下界がその次元に一致するためdenseである。
一次部分のgeneric上界15は別の全関数体rank証明の独立reviewを使用し、今回の標本rank15を上界にしていない。
primitive gaugeへの線形同型により旧weight1座標でのrankとprimitive beta rankは同じである。

同じ実61値からprimitive betaを再構成し、Fisher–Sadek不変量をQQで再計算した。
保存されたc4,c6,Delta,jと一致し、Delta=(c4³−c6²)/1728も確認した。
Delta非零の実点が新しい既約15次元一次像にあるので、像はcanonical-curve discriminant hypersurfaceと一致せず、
一般点のcanonical curveは滑らかである。一次像の別の認識式を推測した意味ではない。

source-recordについて、component ID、全入力名・順序、u1,...,u61の順序、圧縮・展開双方のSHA、
rank16 witnessのpath、次元16/15/1、未構成prime/hypersurfaceのflagsを検査した。
歴史的known51 catalogは変更せず、新recordを加えた集合の52個のIDはすべて相異なる。
既存generic normalizerや設定を再導入せず、受理済みの同じsource hashを再利用している。

fresh producerのmapとsampleもbyte単位で一致した。record内の出力先pathだけは
fresh work directoryからpublic directoryへ対応させて一致を検査した。
独立コードは `code/q_new_family_dense_s61_review_20261103_v63.py`、receiptは
`reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63_REVIEW.json`。
