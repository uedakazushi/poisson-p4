# 全661 sourceにおける未評価honest-second感度

2026-11-05。残る二次adjointのうち、old regular quotientの二次係数だけへの感度を調べる。
これは実高次係数の計算ではなく、既存origin/firstとliteral source ASTの微分である。

真のold r/g leafの二次係数を独立に変えると、pure normalizer由来のm、実C1、X1、hXp0は変化しない。
従って`δX2=Res(A0)δBold2+δP2`、`δrow20 X3=A1 δBold2`である。
row20のRes(A0)は零なので、hXp1とC1にもこのsecond変化は現れない。
18単位blockによる補正は`δa2=Ji δX2_I`、81 blockによる補正は
`−C1 Wi H2(m0 δX2)`を保持する。三wedgeの直接項は
`(5/2)δrow20X3_zir² + mu1 δX2_zir − hXp0 δa2_r²`である。
これらを六dualへ全て合成し、完全661 constant sourceを通してから感度を判定する。

H2は既証のordered-Hadamard/Euler templateそのものである。各sector最後のcomponentのpure-r二次項が零なのは
実cycleのz=0制限から従う。その他の係数や欠けたgauge係数を零に置かない。
新しいimplicit/source/gauge solveはしない。未知感度が残る場合はそのliteral leaf/monomialを保存し、追加係数の値を推測しない。

固定source検査は0.975秒で完了し、さらに強い結果を得た。全1239 typed-B rootsと723 P rootsをsourceへ戻すと、
非零primitive gradientはそれぞれ32 rootsだけであった。それらの二次感度は全661 constant sourceのraw Bold2とP2そのもので消えた。
従ってsix-adjointやH2を使う前に、すべてのhonest-second係数が不要になる。
報告の`all_formal_sensitivities_tested=0`は検査を省略した意味ではなく、raw gradientが消えたためその後のadjointへ渡る非零sourceが無かったことを意味する。
未知のmixed h2/g2を零と仮定したのではない。物理sourceの保存済み二次係数に依存するknown値は別に保持され、まだこのpacketでは評価しない。
