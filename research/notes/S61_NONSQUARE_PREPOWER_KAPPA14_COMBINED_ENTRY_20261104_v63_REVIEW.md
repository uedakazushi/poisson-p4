# κ14と三raw crossの独立QQ監査

2026-11-04。current15の四行[2,4,5,8]、元32の列14だけを対象とし、全73 source演算と
保存16多項式を独立DEGREVLEX算術で確認した。producerはimport／再実行しない。
四κは既受理source、anchor22は既受理bottom24 sourceへ照合し、四M14は実Q2/L prefixの
添字、QQ/a分母、Laurent offsetまで戻す。prefix後の行を裸の135 generatorとは呼ばない。

元の式の符号をsource graphで確認した後、別の逆順四項積でκ14を、各crossは別の積順で

    M[2,22]*M[i,14] − M[i,22]*M[2,14], i=4,5,8

を再構成した。四つ全てQQ多項式として零である。これはκ14だけの零性より強い。

既受理の非零anchor M[2,22]と六固定pure列rank≤1により、分数体上で列14もanchor列に
比例する。従って四行の固定七列[22,11,13,25,15,19,14]の全二minorは分数体で零、
整域への単射から元の係数環でも零となる。anchorを各点で反転する仮定は不要。
このraw4×7 rank≤1は七列に限る。他列、残余rank、全FULLU57・239へ拡張しない。
全追加tailとchart補集合は保持する。新しい点、素数、GB、source jetを計算していない。

再現: `python code/s61_nonsquare_prePOWER_kappa14_combined_entry_review_20261104_v63.py`。
