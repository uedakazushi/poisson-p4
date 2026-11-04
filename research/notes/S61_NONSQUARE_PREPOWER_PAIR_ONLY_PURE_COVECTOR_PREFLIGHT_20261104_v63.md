# 零raw8の後：pair-only四行の直接pure cofactor

2026-11-04。係数・minor・rankを新たに計算しないsource-only packet。
既存raw8は完成して零だった。四行の選択六pure列はrank≤1であり、残る列へは外挿しない。

四sourceはcurrent15 `[2,4,5,8]`、compact14 `[1,3,4,7]`。FULLU既知単元で正規化した
heavy二列をH0,H1、既存pure column22をPとし、A=[H0,H1,P] とする。各rowを除いた
3×3minorの符号付きcofactorをκ_i=(-1)^i det(A_without_i)と置けばκA=0が普遍的に成立する。
他のheavy三列は元sourceでliteral zeroなので、κはactual pure source covectorである。

既存4×6pure blockの全90二次minorが零、Pの非零polynomial entryがあることから、κは
六列 `[22,11,13,25,15,19]` 全てを消す。ここでは非零entryによるdomain cancellationで
多項式恒等式を得るので、そのentryを特殊点で反転せず、追加chartを捨てない。

従ってκ-sourceは残る十一列 `[1,2,3,4,6,7,8,9,14,17,18]`だけにsupportedする。
全8Q9tailを明示的に残す。このsource行ではP4/H5/P5のpivot列が既に零であり、それらの
adjugateやSchur productを計算する必要がない。

Aのgeneric rank3は、既存Hpair rank2とSaux非零から分かる。rank≤2なら、rank一である
六pure列全体がheavy二列のspanに入り、u4のheavy消去後の六pure係数は全て零となり、
既存Sauxが零になるから矛盾する。従ってκはgenericに非零。ただしκの指定一係数を
全FULLU単元とは呼ばない。κ_u7≠0のchartでsource rowu7をκへ置き換えられる。
anchor六列を順に用いるcofactor chartsはD(delta2*Saux)を覆うが、一anchor22だけで
全てを覆うとは主張しない。

このpacketではanchor22に対するκの四成分と、残る十一scalar rootの支持上界、既知B pole
clearing、source atomと型付き演算だけを保存する。符号は全て固定し、値は未評価。
最小のscalar fragmentを後で選ぶ資料であって、full rank/global57や有限bad locusを結論しない。
