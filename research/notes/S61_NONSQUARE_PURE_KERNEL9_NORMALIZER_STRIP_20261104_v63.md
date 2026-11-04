# pure-kernel forward分母の既知単元除去

2026-11-04。固定済みPURE_KERNEL9_CHARTSは変更しない。前向きsource式で逆にする38個の非定数NORMALIZE_gのうち、literal testで未知とされた26出現を取り出す。QQ定数倍を除いて七つ、合計374項である。

一回60秒・2GB・32MiB・150000項の範囲で、a,d,F,B,C,Jだけによる厳密な除法を行う。各段階でbefore=named_unit*quotientを全係数検証し、全積も元のmonic sourceへ戻す。一般factorization、未知因子のGCD、新しいminor・rank・GB・点・素数を使わない。残る多項式はすべてlocalizerとして保持する。

この恒等式はFULLU上のforward recipeの必要分母を小さく記述するものであり、残余因子の零集合をactual cubic defectと同定しない。元の136因子+15pivotのrank同値chartも変更しない。全域の純係数rank・full57・quartic239をこの処理だけで結論しない。

再現: `.venv/bin/python code/s61_nonsquare_pure_kernel9_normalizer_strip_20261104_v63.py --output work/pure-normalizer-strip-fresh`。新しいstrip証明書の独立レビューは別途必要である。

実行は0.100秒で完了した。七入力に対して18回の既知単元除法を行い、各段階と全積を検証した。active係数1928項、全artifactは約350KBである。二入力は既知単元のみとなり、残りは四種類（3・6・3・12項、計24項）のmonic多項式にまとまった。二つの入力は同じed+16d/21+20/3を残す。他の残余にはad−6e/5−151/110、および(ad−6e/5)^2−(226/55)(ad−6e/5)+1104/275がある。これら四つを一般factorizationで分ける計算はしていない。

残る12個のliteral既知単元normalizerと合わせると、FULLU上のforward分母chartは、四つの残余の積が非零という条件で正確に記述できる。どの残余の零集合も取り去ってよいとの主張ではなく、その上でpure cubicを調べる課題は残る。POWER-free17×24から三つのQ sourceを保持して12×19 detectorを作る別の案は、この四分母artifact自体を避けるための次のsource preflightであり、ここでは実行していない。
