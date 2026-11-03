# 零一次ジェットの被約ファイバーと認識の独立監査

2026-11-03。全QQ producerを新しく実行し、保存全係数・認識行列・旧表示照合と一致させた。
別に保存60係数から実際のP=J+Qを組み、全パラメータについてそのJacobi式が恒等零であること、
divergenceと−J grad(trace Q)の一致を直接再計算した。
固定線の全60係数認識行列のrankは50,50,51、像基底は10,10,9であること、
積が零であること、全13/13/12次元の微分minorもQQで確認した。
旧v07の九パラメータ表示とc=0,1/4のラベル、およびR0のconstant-factor定義を読んで照合した。
R0の次元12と、それを含むM0/M1の次元13を混同していない。

全CA/CB/Nからのproper全射をβ=0へ制限すると、任意のtarget点は零基底変換へ持ち上がる。
六つのβパラメータは左逆によりschemeとして零になり、sourceはrank10/10/9束である。
従ってproper像は全underlying reduced zero fiberであり、密な族だけの主張ではない。
この段階でtargetのscheme-theoretic zero fiberが被約とは推論しない。

δQ≠0ではそのrank1像が旗線を一意に戻し、旗方程式の線形条件が厚みも消す。
δQ=0では固定線の全認識式がv∧Q=0だけとなり、旗ファイバーはschemeとしてP F(Q)である。
非零Qの定数因子空間は一つの非零decomposable Q(y0)の二平面に入り、次元は高々2。
従って被約一点、P1、原点でP3という全source旗分類を得る。
両Mの非零δ点では線が共通なので、固定線の交叉からM0∩M1の被約台はR0となる。

内在的な有限式のirrelevant saturationと消去によるprime表示も正当である。
全固定線のkernelが定数rankで、G同変性により射影入射はintegralな線形束となる。
vのirrelevant idealで飽和したconeのidealはintegral入射の閉包のprimeであり、
そのq係数環へのcontractionがprimeとなる。proper性が像を閉にする。
展開済み最小生成系やtargetの非被約primary idealを得たという意味ではない。

rank Kβ=2の診断は、真の基底入射には追加収縮条件があるためP ker Kβではないことと整合する。
この点は全R1監査の普遍的旗一意性論証でも既に確認した。
全S61の零ファイバーや全target scheme構造の尽くしは残る。

検算は `code/ambient_zero_firstjet_review_20261103_v63.py`、
結果は `reports/AMBIENT_ZERO_FIRSTJET_FIBERS_20261103_v63_REVIEW.json`。
