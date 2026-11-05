# 同じ列17の十二角軸入力：追加sourceだけのprofile

2026-11-05。新しい軸係数・determinantは未計算。三Cramer独立peerによりrawpairの
列17比例が確定したので、同じr17=det(M12)/P_K92というsource恒等式を使う。
lambda多項式テストは不一致のまま保存し、候補十二項を実lambdaには使わない。

既存M11の軸Caxis全579整数係数を再利用する。generic M11に行一＋二行二という
determinant一の前処理を施したのが旧Mprimeなので、追加する列17にも同じ行操作を施す。
新しいG134行は前処理しない。列順は[22,11,13,25,15,19,14,18,28,30,29,17]。
最後の三H3列はL=aFC²で割る。G134のrho13=(8000000/16119)Lなので、追加H3三値の
sourceはcompact13の列[1,3,2]の定数倍であり、F/Cの逆多項式を新たに作らない。

追加pure行と列は元Q2/L finite DAG、Q9は既存h9 literal sourceを使う。
checkpoint全体のJSON objectをloadせず、必要なmatrix25x32 cellだけをdecodeする。
sourceの分母はQQ*a^k、またはQQ*B。B(a,0,0)=-40で非零なので、各追加行/列を
既知単元aとBの冪で払えばe,d非負の多項式となり、軸制限は合法である。
軸d0をFULLUのPoisson点とは扱わない。kappa3や未知chart多項式の逆元は入力に入らない。

source atomのe=d=0支持だけを数え、L1・QQ分母・a区間を和積の上界として伝播する。
signed係数列は生成しない。既存row/column a-shiftを保ち、新row/columnだけに必要な
追加a-shift、B-clear、row LCDを選ぶ。旧十一角を含む全十二角整数化のdeterminant係数と
次数・高さ・owned・RAWの上界を保存し、実workerのGOとは区別する。
構造上のaxis zeroはgeneric zeroを意味せず、非空の支持上界も非零証明ではない。
