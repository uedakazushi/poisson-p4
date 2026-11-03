# V9 の二つの短い単元 pivot の独立監査

2026-11-03。`S61_BINARY_NONSQUARE_SHORT_UNIT_PIVOTS_20261103_v62` の二つの定数列操作と
宣言された単元開での一段の rank 分離を承認する。全点 rank41/57、cokernel の自由性、
quartic65、ambient の尽くしは依然として未証明である。

凍結 archive の全 member、STATUS、既存 V6・PRODUCT_ALT への依存 hash を確認した。
さらに現存する元の全14×46 product 行列と約89MBの V6 P-clean27 を stream で hash し、
その元データの row10・row4 を読み、archive の target row と全列で一致させた。
target row が別の生成段階の行列から取り違えられたものではないことを確認している。

探索コードや native polynomial 演算を呼ばず、各疎多項式を exponent→Fraction の辞書として読み、
二つの保存 witness の全61/46成分を直接掛け合わせた。その結果は正確に

    (−14985/8)M[10,11]+(2997/8)M[10,13]=a³(5ad−6e−4)²,
    −(N[4,5]+N[4,30])/50=(100a+18e²−9e−2)(5P−6e−4)

である。main の三番目の serialization 変数 d は P=ad を意味する。
保存された六つの単元を実際に読み、その積の指数から同じ二つの多項式を再構成した。
新たな因子を単元として加えていない。

同じ列の一次結合をすべての行で実行し、列11または列5を置き換える行列の行列式は
それぞれ −14985/8、−1/50 である。どちらも非零有理数なので、列操作自体は QQ 上で可逆である。
その target entry は宣言された開上の単元であるため、通常の行・列消去で1×1単元 block と
残余13×45または26×60へ分かれる。native 実装が行全体から除くものも、非零有理数と
明記した六単元の共通因子だけであり、この局所化上では可逆操作である。
従ってすべての点で、元の rank は1に残余 rank を足した値に等しい。

約94MBの main 残余をこの監査で再生成したとはしない。残余 term 数44143/1290701と SHA は
凍結 native 実行の receipt を使用し、今回の独立計算は元の全 target row、完全な短い恒等式、
単元因子、列操作の可逆性に限定する。新しい一段を既存の単元消去と合わせると、
main は31段、mod-ell は11段の構成的な単元 pivot を持つ。
全 residual rank や点別 cokernel をこれから自動的に決めることはできない。

独立コードは `code/s61_nonsquare_short_unit_pivots_review_20261103_v62.py`、
receipt は `reports/S61_BINARY_NONSQUARE_SHORT_UNIT_PIVOTS_20261103_v62_REVIEW.json`。
