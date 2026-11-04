# 正weight Bianchi のweight4包含：独立加群監査

2026-11-04。対象は frozen
`S61_NONSQUARE_BIANCHI_WEIGHT4_CONTAINMENT_20261104_v63`。
新しい係数展開・rank・pivot・元175恒等式の検査は行わず、受理済み source
と whole-module-kernel の二つの定理を結ぶ。

`R0=QQ[a,e,d,1/(aFBC)]` は整域である。元63 cubic sourceの六つの厳密関係は、
削除順48,55,59,57,61,14に対し、pivotが1,1,1,a,−a,12a/5。
後の関係に既削除行がないので、その単元三角 minor は −12a³/5 である。
従って六関係による source 商は自由57。元の91列への写像の generic rank57
により、この自由商の kernel は generic に零である。整域上の自由加群の
部分加群は torsion-free なので、この kernel は零。ゆえに六関係は **R0 上の
kernel 全体**の基底となる。特殊点での周囲 rank57 や cokernel flatness を
用いていない。これは既受理 row-basis レビューと同じ論理である。

full operator の保存された全 weight pieces に3未満はなく、weight3 piece
は indices35,56,71,81,87,90の六つだけ。いずれも spatial degree5から来る。
選択34式の exact graph と恒等零 nuisance19/23を除くと、6 quadrics×weight1
coordinate と21 genuine cubics×constantだけに支持を持つ。weight1 graph は
T·D·vで尽くされる。従って既受理 full135 identity と frame/row-unit 共変性に
より、六 piece はすべて同じ63 cubic sourceの R0-kernel の元となる。
それらを六基底で表す係数を新しく求めなくても、whole-kernel 定理が包含を保証する。

weight4へ到達する正weight polynomial multipliers は、weight4 pieceへの
constant、weight3 pieceへの七つのlightのlinear combinationである。
H はweight2なので後者には寄与しない。六基底の最初五本×七lightは、
独立受理された全35単元削除と一致し、354→319商で消える。第六基底は
元0..41 product rowsだけに支持を持つ。そのlight倍は六quadのlight²-cofactor
関係であり、H係数は現れない。従って旧168-column純quadratic mapの kernel
に入る。この **whole module kernel が旧22の自由加群に等しい**ことは FULLU
上で受理済みで、特殊fiberのrank146を仮定する必要はない。

E22を旧22、D15を新15 direct weight4 pieces、B4を今回の operator pieces
から通常の正weight多項式倍で作った weight4 source関係の319行への像とする。
上の議論により

    D15 ⊂ B4 ⊂ E22+D15,  E22+B4=E22+D15.

既受理の I11/I4 は D15 の15本を aFBC 上で自由かつ点別独立にする。
E22は quadratic source列だけに支持を持ち、この15列では零なので交わりは零。
従って **FULLU 上で E22+B4=E22⊕D15≅R^37**。80方向 gate には43本不足する。

これは **B4単独のrank37ではない**。第六cubic関係のlight倍が旧22のどの
部分を張るかは今回求めていない。higher-weight pieces、monomial除算、
飽和、負weight操作、別のsyzygies、full80、C122 rank64、全点quartic239、
affine exhaustionをこの包含へ含めない。tG6は受理済みのままweight6にあり、
tを単元として落とすこともしていない。

専用 checker はこれらの受理済み hash、六pivot・source support、35削除の
labels、15個のdirect piece、二つの whole-kernel gate を固定する。新しい
数値点・field・GB・大行列の算術は含まない。
