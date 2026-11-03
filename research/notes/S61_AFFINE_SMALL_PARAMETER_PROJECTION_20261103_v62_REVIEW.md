# 9パラメータ階数・8パラメータunit判定への削減の独立監査

2026-11-03。S61からprimitive二形式への線形写像、24×6断面の一致、
16次dominance determinant438994543750/3をQQで独立に確認した。
元134式の重み70本の4、64本の5も照合した。
さらにCASの多項式環を使わないFraction疎多項式積で、五変数断面への
全134式の記号的代入係数を別に構成し、保存全式と一致した。

Sp4×Gm×A5のdominanceとq上のaffine線形作用により、これは一般βの特殊化ではなく、
支配的体拡大とfiber同型による正確な削減である。
二つのmoment-curve線形式の中心回避は、真の射影境界の次元≤1と、
各非零qに対する二多項式の有限零集合、proper射影の閉性から従う。
基点のない[t:Lz:Lw]はO(1)のamplenessにより有限で、全二次元cycleの次数がgeneric module rankとなる。
次数30の定理は元のC(W16)-整曲面に適用する。
体拡大でその曲面が分裂しても総cycle次数は保たれるため、rank≤29なら曲面を除ける。
断面体上の個々の分裂成分へ次数30を課してはいけない。

Wronski対は別のunit判定である。
曲面関数体の二つの独立な微分から得る係数多項式P,Qは独立で、
標数零ではPQ'−P'Q≠0となる。これはLz,L'zのJacobianなので両者は代数的独立。
従ってgeneric二値fiberがunitであることは、次元≥2のsupportがないことと同値である。
この対にprojective degreeとの等式やrank29判定を適用しない。

二線形式のq0,q1係数minorはsecantでw−z、Wronskiで1であり、実効43未知数を保つ。
九／八パラメータのsolverのrankまたはunit結果は今回得ていない。
曲面除去、曲線除去、五点尽くしの証明とは区別する。
検算は `code/s61_small_parameter_review_20261103_v62.py`、
結果は `reports/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62_REVIEW.json`。
