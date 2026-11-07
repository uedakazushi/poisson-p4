# 既存quadratic_A λ7点は滑らかで唯一のMPS24成分を通る

2026年11月7日。producer exact squeeze、独立レビュー待ち。
全七源のMPS帰属とは別の、同じ既存点における強い帰結である。

点は旧 `CLASSIFICATION_NIL22_QUADRATIC_A_WEIGHT_FRAME_20261107.json` と同一の
λ=7、m=(1,0,0,0,0,0)、整数代表2Pである。
元の一般witnessのI=8を正規化した同じ軌道であり、新しい点を探索していない。
旧全350×150微分Dのrank101、Euler rank24、実全七源を含むF rank48を入力とする。
従って射影Poisson schemeの全接空間はdimension24である。

`QUADRATIC_A_SYMPLECTIC_MPS_LIMIT.md` のliteral ε曲線を元座標へ戻し、
同じ整数倍2を取り、そのε一次微分Rを追加した。
producerは全350の `D R=0` を直接再確認した。
全25Euler、全25GL、scalar、**全七source**微分とRを並べた150×59行列は、
QQの非零49×49minorによりrank49となる。
旧Fは同時に再構成してrank48、旧点も全150係数で一致することを確認した。

このrankは、無関係な二族の接空間を足しただけの下界ではない。
全λ ε曲線に、全六磁気normalizerの逆PGL変換とinverse tensor scaleを作用させる。
αβIはこの同じ点でunitなので、これは点の近傍で正則な一つの実joint族である。
ε=0では元の全七sourceそのもの、m・λ微分は旧七source方向、ε微分はRとなる。
さらにPGL5作用とscale/Euler自由度を加えた実像微分が上の行列である。
従って標数零の実像次元下界はrawで49、射影Euler商後で24である。

全接上界24と実像下界24が一致するから、元Poisson schemeはこの点で正則であり、
通る全空間成分は唯一、dimension24である。
独立受理済みの全源帰属曲線により、その成分はMPS40の一つである。
特定のT番号はここで同定しない。

旧十八変数模型は真の48本cutと完全IFTによる局所同型であった。
新actual方向が旧F48の外にあり、全局所次元がraw49なので、そのcutの被約支持はdimension1である。
したがって、以前停止した「normalの純冪をidealに入れて支持を原点へ縮める」という目標は偽だった。
旧Mora85秒・unit90秒の未完了結果に、過去に得なかった結論を書き足すものではない。
新しい実曲線と次元一致によって、現在の数学的状態を訂正したものである。

将来、別の未被覆源の既約実族がこの点を通ることを**実際に**証明できれば、
成分一意性によりその族全体を同じMPS成分へ運べる。
単にその族の境界因子型が一致することからは適用しない。
double_AについてD3を零にする素朴な固定N源はV1のtraceが−3/4であり、
quadratic_Aのtrace−1に一致しない。未証明のλ=1/4橋は採用せず、追加point検査もしていない。

exact producerは `classification_quadratic_a_smooth_anchor.py`、
reportは `CLASSIFICATION_QUADRATIC_A_SMOOTH_ANCHOR_20261107.json`、
jobは `jobs/quadratic_a_smooth_anchor.json`。全行列・minor・R・入力hash・実行bytesを公開した。
計算0.608秒。最初の読込は保存parameter名lambdaがPython予約語であるため停止したので、
parserで同じsymbolのalias lamへ直して再gateした。数学・点・sourceは変えていない。
