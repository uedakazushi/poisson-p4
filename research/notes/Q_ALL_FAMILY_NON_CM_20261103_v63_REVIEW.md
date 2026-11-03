# 全 Q_a の full Poisson 完成局所環の非 CM 性の独立監査

2026-11-03。新規 producer report は保存結果と完全一致した。
proper CA/CB 模型の raw tensor を補助 builder なしに明示式から再構築し、
両者の一致、divergence 零、座標順 (old0,old1,old4,old2,old3) による
physical Qa の全行列成分の一致を QQ[a] 上で確認した。
CB の比 [a:−2/9] は全有限 a において定義され、Qa は a=0 でも零でない。
従って閉包や rank4 開の特殊化を推測せず全 a で CA と CB に入る。
実 T01 帰属は別の独立監査済み可換行列写像による。

ここで用いる25次元 T01、28次元 CA/CB の full Poisson 空間における
実際の既約成分性は、既存の component 証明書の結論である。
異なる大域 minimal 成分は閉点で局所化しても異なる minimal 成分として残る。
有限型 C 代数の優秀性と完成の忠実平坦性により、各成分の25または28という
局所次元を持つ completed minimal 成分が存在する。
完全 Noether 局所 Cohen–Macaulay 環の unmixedness に反するため、
全有限 a で full Poisson 完成局所環は非 CM である。
同じ異なる minimal 成分は被約化でも残るので、被約空間を完成した環も非 CM である。
nilpotent の存在はこの証明の前提ではない。

個々の CA/CB の全 a における非 CM 性、全局所枝、局所次元の上界は結論しない。
監査結果: `reports/Q_ALL_FAMILY_NON_CM_20261103_v63_REVIEW.json`。
