# 実heavy minorだけを使うpure12×19入力案

2026-11-04。保存済みPOWER-free17×24は三つのQ source（旧label4,5,9）、14genuine、19pure target、5heavy targetである。pure323成分のparameter項は70399（Q3は692項）、heavy14×5は既受理の70個のFULLU単元scalingで2407項のcompact Nへ移る。合計72806項で、旧pure9 chartの140820項より小さい。今回はsource/profile/型の照合だけで、新しい係数・minor・rankを計算しない。

actual H7上でNのrank5を使い、五つのgenuine行Iの**実際の**heavy minor Delta=det(N_I)が非零なchartを選ぶ。残る九つのgenuine sourceをDelta*e_u-(N_u adj(N_I))*e_Iにし、三Q sourceをそのまま加える。この12 sourceのpure19像が12×19 detectorである。source basis minorはDelta^9で、rank(17×24)=5+rank(12×19)。rank12はfull cubic57の十分条件だが、その不足はfull57の否定ではない。

新しい分母はDeltaと既知FULLU単元だけであり、degree207 POWER normalizationも旧四つのNORMALIZE残余も導入しない。五heavy行を選ぶ全candidateのうち、非零項を持つ11行から5行を選ぶ462組をsupportだけで調べると、371組が構造的には非零minorを持ち得る。これらのdeterminantは計算しておらず、全て非零だとは主張しない。rank5ならこの有限候補のどれかが非零なので、chart coverとしては十分である。

固定sparse pairはcompact行[1,3]、残りのB7行は[6,8,9,10,11,12,13]である。D(p)ではpairdet=-(6/5)*B*p²が単元なので、Bの三行minor35個が五heavy chartを与える。ただしp=0でもH7が成立する部分があり得る。BgoodをFULLU全体でH7と同一視せず、このsurfaceや別のactual H minorを保持する。

## さらに小さい既存sourceからの代案

既受理の二つの**global pure Q** completionだけを使うと、POWER前の15×22（Q9一行+14genuine、17pure+5heavy）がある。ここで五heavyを消すと、Q9を残した10×17になる。Q9の九つのactual係数h_iは既に373項・次数(6,7,2)・最大41bitで明示されている。既受理のsum H_i h_i=a^29とa単元性は、九つのD(h_i)がFULLUを覆うことを与える。高次cofactorを再展開する必要はない。

各D(h_i)でQ9を一つのpure targetへpivotすれば9×16になる。この際、九つ以外の八つのpure tailも全て保持する。必要な分母はactual heavy minorとactual h_iだけで、旧normalizerの四つの人工的chart因子を避けられる。rank(17×24)=2+rank(15×22)=8+rank(9×16)である。ただしこれは同じprojected detectorのrank同値であり、元のfull cubic cokernelとの同値だとは言わない。

新しいminorや係数を作るには別の有限予算が必要である。再現: `.venv/bin/python code/s61_nonsquare_pure_kernel12_source_preflight_20261104_v63.py --output work/pure12-preflight-fresh`。新点・素数・GBは使わない。

source preflightは0.103秒で完了した。全70 heavy entryの旧checkpoint一致、三Q行のheavy零、compact B7×3への全21 entry一致を確認した。371候補の列挙はnonzero supportのmatchingだけで、parameter値もdeterminantも計算していない。373項のh9 sourceを新しい小artifactへ抜き出したため、後続のsource構成は高次POWER cofactorを評価せずに使える。actual h9とPOWER前の15×22行列の対応は既受理のPOWER source reviewを明示依存として保持する。
