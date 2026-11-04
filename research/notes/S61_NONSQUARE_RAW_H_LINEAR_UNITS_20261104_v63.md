# 元の16三次式のH一次係数を使う15×6 criterion

入力は凍結済み57sourceの最後16 genuine cubicのH*light係数112個。大きなH1
quadric-cofactor行列を再評価せず、この16×7だけを使う。physical lightから
(x0,x1,x2,x3,x4,Z,W)への定数frameは既存のexact7×7変換と逆を確認する。
源C6をCdelta=delta*(C6−C1)、delta=−64000C²へ置換する。

そのH係数はdelta*z Z+delta*w Wのみ。既存のliteral Bezout
cZ*delta*z+cW*delta*w=eta=192aFB²C³/125がある。従って係数columnのGL2
`[[cZ/eta,cW/eta],[-delta*w,delta*z]]`でCdeltaはunit座標一つになる。
第二行を既知単元gamma=delta*F*B/10^7で割ると、残る15×6は多項式で
五つの元係数と`h*ci_Z−a*g1*g2*ci_W`からなる。etaを展開して全15行の
分母へ掛けない。捨てるunit座標の係数はnumerator_i/etaとしてlazyに保持し、
同じ実際のfull cubic sourceからCdeltaを引く。pure部分を捨てない。

支持が10000項以下の場合だけ、既知a,d,F,B,C,J単元を持つ個々のentryをpivot
とする一回のGaussianを行う。QQ定数と全行共通の既知因子によるrow清掃だけを
許し、未知因子を割らない。各完成pivotで操作と残行列をcheckpointする。
元行列の選択minorの単元式は各rowの累積単元scaleを追跡して復元する。

もし残6列全てが単元pivotで消えれば、元H-linear16行はFULLU全域で7次元を
生成する。そのとき全112のlight倍はH*Sym²(light)28方向を生成し、既存target
へのH1射影23も全射になる。full sourceのpure係数は同じrow操作でlazyに残る。
この含意をfull57/full239へ読み替えず、source自由生成元数の自動削減も行わない。
止まった場合は得たunit lowerだけで、全点rank defectや非全射を主張しない。
