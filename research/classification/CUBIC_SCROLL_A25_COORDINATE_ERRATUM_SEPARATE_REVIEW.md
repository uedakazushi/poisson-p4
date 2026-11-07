# A25 coordinate erratum の別作者確認

2026-11-07、reviewer `/root/log_manuscript`。純手計算 **PASS**。
対象は `CUBIC_SCROLL_A25_COORDINATE_AND_POLE_ERRATUM.md`。
producer、旧review、旧manifestのbytesは変更しない。新CAS・source計算は実行していない。

a=x2=1では二incidenceは

`x3−b+c x0=0`, `x4−b x3+c x1=0`。

第一式を解けば `Y1=x3=b−c x0`。
これを第二式へ入れると `Z=x4=b²−bc x0−c x1`。
従って

`Q=Z−Y1²`
` =b²−bc x0−c x1−(b²−2bc x0+c²x0²)`
` =c(b x0−x1−c x0²)`。

両incidenceとquadric equationの全literal恒等式を確認した。
generic p^*LではQのc-orderは一、c=0でのcollapsed conic coordinateはY1=b。
この二つが元の全pole-order/residue proofで用いた性質なので、正しいminus符号へ直しても結論は変わらない。
full fibers、重複kappa、c|Omega全補集合の扱いにも変更は無い。

もう一つのpole-location訂正は、既凍結 `CUBIC_SCROLL_A25_BOUNDARY_SEPARATE_REVIEW.md` §6で独立に確認済み。
reducible Omegaではcomponent交点のpoleを許し、smooth irreducible Omegaを得た後だけ
poles⊂L∩Omegaを使用する。
上の正しいcoordinate式とこのscope補正を合わせ、smooth cubic scroll rank二の全G像≤24という
producer結論と既boundary-reviewのPASSを維持する。
