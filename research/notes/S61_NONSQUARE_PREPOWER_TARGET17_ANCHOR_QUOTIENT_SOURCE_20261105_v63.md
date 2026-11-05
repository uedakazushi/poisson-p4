# 実 anchor 比と同じ列17の十二角 Schur 分子

2026-11-05。三Cramer primary [0,0,Y2] の独立peerは別gateとして保持する。
κ17零は独立確認済み。Y2は1768項の非零多項式であり、実比例係数は

    lambda = -Y2/kappa3.

κ3は既存chartの単元。literal G20,22成分を新たな単元とは主張しない。
旧anchor列のpivotが既に除去されていることと、ある一成分の単元性は区別する。
全四pair行について col17=lambda*col22 が従うため、列17を
col17-lambda col22に置き換えると、その四成分は全て零になる。
十角補正は省略しない。既存四＋十の同じ列は

    corrected_v10 = v10-lambda*C[:,anchor],
    corrected_v1  = v1-lambda*c[anchor],
    r17 = corrected_v1 - corrected_row10*S10inverse*corrected_v10

となる。全H3、source正規化、十角inverseとtail2の既消費効果を保持する。

さらに新しい係数計算なしに、同じr17を十二角borderで表せる。
既存十一角M11はrow順Q9,G18/a,G39/a,G68,G20,A,B,D,C,E,Fであり、
列順[22,11,13,25,15,19,14,18,28,30,29]。最後の三H3列はL=aFC²で割る。
これにcurrent15の行14=G134を最後の行、original32列17を最後の列として加える。
新しい十二角行列をM12とすると、既存chart上で厳密に

    r17 = det(M12)/P_K92.

以前は不正だった裸の十二角borderを、この式では強いcol17比例の証明後に使う。
既証八pure列と新しい列17の全てでraw pair四行のrankは一以下、三H3列では零。
generic fieldで一時的に非零anchor成分を反転し、他の三pair行をG20の倍数で消すと、
selected十四角・十五角の両determinantに同じ三角heavy/tail2因子が現れて相殺する。
三行と三列の同じ移動の符号も相殺する。H3のL三乗、G18/G39のa二乗、
A/B/D/C/E/FのQQ単位とa正規化も両determinantで同じなので相殺する。
G134行と列17はscale1。結果の両辺は元のchart上で正則なのでdomainへの単射から
恒等式は下降し、新しいanchor逆元や補集合の削除は不要である。

det(M12)は未計算で、r17の零/非零も未決。他の七列、全補集合、全FULLU57は保持する。

lambdaを更に小さくするSOURCE-only可能性を固定する。もしlambdaが多項式なら、
各変数のmax/min次数と全次数の加法性から、その支持はa0..2,e0..2,d0..1、
全次数四以下の17単項式に入る。これを仮定して結論せず、17個のscalarだけを
固定leading termの三角式で解き、全係数でkappa3*lambda+Y2=0を確認する手順を
次の候補とする。失敗なら既存の有理式を残し、別candidateへ自動拡大しない。
多項式GCD・factor・新点・新素数・別targetは使わない。実17係数は未計算である。
