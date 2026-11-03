# 実交叉の全有限表示と平方零 nilradical：独立監査

2026年11月3日。対象は S61 内の実 target intersection D=CA∩CB の指定 a=1 Q における完成環。
共鳴 a=2/9 の全35 model、全 P4 の rank≤2 境界は対象ではない。

frozen builder を新しい作業 directory で QQ 再実行し、秒数以外の全証明書が一致した。
最終 M2 colon/reduction script も独立に再実行した。
各109式が全 kernel となる理由も確認した。定数可逆な target K block から source fiber 全体を
多項式で戻し、偶奇分解後は T⊕Jp の像となる。J=(a,s,7 independent linear forms)は regular sequence。
pair products による全式の次数低下と、J の全第一 syzygies が Koszul で尽くされることから、
28 linear+36 syzygy+45 product=109式は全 image kernel を与える。
必要条件の pullback 零だけから kernel equality を推測していない。

追加の独立コードは、N の実全61座標と九つの実 target K を使い、N の全 kernel を別に構成した。
こちらは JN=(a,s,6 independent linear forms)で、30 linear+28 syzygy+36 product=94式である。
全94式の source pullback と regular-sequence rank8を確認し、新しい M2 script で

    ID+(F)=independent_IN=IA+(F,G),
    ID:F=(all49 coordinates except v1,v4,v5,v7)

をそれぞれ両方向へ検査した。F≠0,F²=0,G=0,v1F≠0,v1²F≠0 も fresh GB で確認する。
従って nilradical は正確に(F)、冪零次数2、全加群は smooth11次元支持上の C[[u11]]F である。
10次元共通非正規層だけが支持だという初期の候補は偽であり、この監査もそれを用いない。

IN は独立に示した integral image の prime ideal である。CA source の F は文字通り rA であり、
(F,N_K9) の全 fiber Jacobian determinant は−72である。
全像 A⊃R′,p²B における A∩rAB=(rA,rAp) という restriction の議論とも一致する。

正確な完全列 0→R/Qemb→R/ID→R/IN→0 から associated primes はこの二つに限られる。
Qemb は実部分加群の annihilator なので associated、IN は唯一の minimal prime なので associated。
完成後も IN は受理済み A7 の domain、Qemb は regular domain であり、他の associated prime は生じない。
七群方向を戻すと左項は正則11次元、右項 A7 は depth11である。
local cohomology の H^10(A7)=0 と H^11(left)≠0 から depth D=11で、dim D=18となる。

218式による全完成環表示、主平方零 nilradical、全 annihilator、正確な二 associated primes を受理する。
環の trivial extension への分裂は結論しない。全域の rank≤2 境界の交叉も依然この局所結果の外である。
証拠は code/ambient_intersection_exact_review_20261103_v63.py と
reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63_REVIEW/ に保存する。
