# lambda固定17係数：共通分母による二千bit以内の源上界

2026-11-05。旧 SOURCE の3524bitは二つの大きい未知係数を掛けるとした粗い上界だった。
本方式ではその積は存在しない。原 SOURCE は変更しない。

原κ3とY2を整数primitive K,HとQQ単元に分けるとlambda=-(H/K)/3。
leading係数の絶対値をb、K/Hの整数L1をk,hとし、D=b^17を全17係数の共通分母にする。
候補係数の整数分子R_mは固定graded lex順で

    R_m=(D H_(m+LMK)-sum(previous_n R_n K_(m+LMK-n)))/LC_K

により求められる。Dの選び方から全17除算が整数exactになる。
これは多項式lambdaの存在を仮定せず、固定候補を作る有限scalar計算である。

A=h(b+k+1)^17と置く。第j段の既計算分子のL1は
h b^(17-j)(b+k+1)^j以下、従ってA以下である。
各scalar除算前の分子、全検証積K*Rのpartial係数、D*Hとの差の絶対値は
D h+17 k A以下になる。source整数からの実bit値はDが998、Aが1651、最後の上界が1746。
R同士を掛けないので倍の3524bitは不要で、全明示整数temporaryは1746bitで足りる。
実装には余裕をもって2048bitのbefore-commit guardを置き、越えればprefixを保存して停止する。

全係数検証の場所は全次数34以下の7770点に限る。17個の三角scalarと、K1174項×17、
H1768項を全て使い、D H-K Rが全係数零の時だけlambda多項式と認める。
失敗ならcountermonomialと候補を保存し、既存有理式lambda=-Y2/κ3を保持する。
候補拡大、多項式GCD/factor、点、素数、別targetは使わない。

packed keyと係数を別slotとして、入力5884、検証15540、17scalarの両バッファ34、
metadata15000と256余裕を含め36714。原JSONを解放してから検証mapを作る。
canonical QQはgcdによるexact整数除算で既約化し、巨大な分母同士のcross積を避ける。
guard prefixは全有効値とcursorをbinaryに保存し、入力hashとparserを固定する。
RAW上限は旧6.82MBを維持する。実17係数・最終一対八はまだ計算していない。
