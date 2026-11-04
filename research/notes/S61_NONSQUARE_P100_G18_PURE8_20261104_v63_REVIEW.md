# P100 の full G18 source・pure8 fixed minor の独立監査

G18_PURE7、PURE8_RANK_SOURCE_PREFLIGHT、PURE8_GENERIC_RANK を監査する。full saved
57×91 の compressed matrix を hash まで確認し、basis position41（old normalized63 row42、
actual135 row18）の全91係数を読む。七つの H係数の literal 零性と、全84 pure core係数を
直接比較する。旧 restricted ZW functional を full cubic として持ち上げない。保存行は既に
Tn normalization 189a/F で割られており、row unit を再び掛けたり割ったりしない。

既に独立受理された peer の QQ[z] named-unit arithmetic だけを再利用し、P100全曲線の
84係数、physical/joint inverse49と両方向の積、joint84、七つの actual source products の
全1470 quartic coefficient identities と全686 NF98 outputs を別再生する。全122 NF output
の H-degree dependency も再構成し、pure source に H1/H2 component が入らないことを確認する。
新 oracle query・source point・rank search・GBは行わない。

q5 cleared row は先の independent pure98 source receipt を保持する。全八行の既知五因子
LCD を確認し、指定 rows[0,2,3,5,7,6,1,4]、cols[0,1,2,3,4,5,6,9] の64成分を literal source
へ結び戻す。作者の border search を繰り返さず、別 FLINT integer polynomial determinant
algorithm で固定8×8の determinant を確認し、保存 adjugate の両側128積も全係数検査する。

非零 delta の上では actual pure8 が独立なので、既存117 quadratic・23 H1・1 H2と合わせ
元 weighted quartic の下界は149。delta の未知零点は除いておらず、全P100下界は142のまま。
q5行の受理済み単元 dual と explicit GL2 completion により original8 rankは1+residual7×97 rank
であるため、残余generic rank7も従う。22MBの残余行列を展開・再生しない。uniform149、full239、
他の曲線や全FULLUへの結論はない。証拠と入力hashはreportに保存する。
