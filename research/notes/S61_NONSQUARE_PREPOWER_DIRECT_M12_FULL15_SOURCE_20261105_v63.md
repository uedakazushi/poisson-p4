# 十二角borderからP chartを経由しない十五行の十分条件

2026-11-05。三Cramer独立peer20f835により列17も含め九pure列でrawpair四行のrankは
一以下である。行列sourceとH3 carryはM12_AXIS_SOURCEとANCHOR_QUOTIENT_SOURCEに固定した。
scalar determinantは計算しない。

Aをliteral current15 row2・original32 column22成分とする。Aは非零多項式だが、
FULLUで単元とは言わない。M12のsource十二行に含まれないpair三行current[4,5,8]を
G134の前に挿入し、M12の最後の列17の前へ[27,31,2]を挿入した十五角をD15_normとする。
全source row単位と列H3/L、H2/L,Lgamma、tail2/Bのcarryを既存A4と同じにする。

すると同じ列・行の全sourceに関し

    A * D15_norm = det(A4) * det(M12) = -T * det(M12)

が成り立つ。一時的なD(A)上で他のpair三行をG20の倍数で消すと、それらは
三つのheavy/tail2列にだけ残る。十四角/十五角のblock Laplaceの符号は三行移動と
三列移動の両方で相殺する。A4でanchor22列を最初に動かすのは二回の交換で符号正。
分母を払う同じ証明では、pair三行をA倍してanchor係数倍G20を引くことで
A三乗が両方に現れ、共通A二乗をdomainで消して上式を得る。

ここでD15_normは元十五行・選択十五列minorのknownFULLU単元倍である。
source orderは[0,1,3,6,2,11,12,9,7,13,10,4,5,8,14]と各固定QQ/a row単位、
target orderは[22,11,13,25,15,19,14,18,28,30,29,27,31,2,17]。
G18/a,G39/aと行一＋二行二の前処理、E=aG129/c12などをM11と全く同じにする。
新しいκ3、anchor、P_K92 inverseをこのidentityに入れない。

従って旧Omega上でdet(M12)を単元とすれば右辺が単元となり、AとD15_normも単元。
元prePOWER15×22は全行rank十五で、既存full57 source liftの十分条件を満たす。
新しい十分chartは旧Omega∩D(detM12)でよく、余分なD(P_K92)は必要ない。
旧Omega因子の更なる削除は別source監査を要するので本記録では行わない。
detM12、全残余、補集合、全FULLU57/239は未計算・未解決のままである。
