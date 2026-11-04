# NEW25 六変数 source・二共変量回路・一意な tangent kernel の独立検証

2026-11-03。著者 producer と native worker を再実行しない。既受理の元の
全三変数 Poisson tensor から n1=1、評価点 `(X,0,Z,W)` の値と全一次係数を
別に作り、全 raw30 成分を照合する。既受理の有限 parameter gauge は tensor
scalar を伴い、modular Ga による y=0 gauge の点移動式も直接照合する。

rational Darboux frame、全24 primitive 係数、受理済み16基底への復元、各
scalar c 指数を有理恒等式として再構成する。保存された1025-node 回路の
beta16 祖先を QQ(u,v,c,X,W,B) 上で評価し、C,D の全二次係数と行列積・trace・
Newton の有限 word により F4、Inv8、pencil5、Delta8、N12、m を束縛する。
pencil5 と trace8 の変換定数は、元の全 W16 における分母を払った整数係数
多項式演算でも独立に確認する。N12 や m の source 全面展開や標本照合は
使わない。beta の open は cB≠0、元の一般族 chart には uv(32u+v)≠0 があり、
m には Delta8≠0 を追加する。分母を恒等単元とは扱わない。

kernel の30×26行列は raw tensor の添字作用で新たに作る。六 source 方向は
`du,dv,dX,dW,dB,dlog(c)`、残りは GL4 の16方向と perspective4方向である。
c≠0で全 raw 行を c で割るため、行列は QQ(u,v,X,W,B) 上となる。
全20個の有理数単位 pivot を照合し、残る10×6行列について**別の明示的な
5×5非零多項式 minor**を計算する。全26成分の保存 kernel は元の30式で
恒等零と直接確認する。従って階数25と一意な一般 kernel は native rank
receipt に頼らず従う。source 部分は非零で du=0 である。

本稿は source 六変数の有理回路と tangent kernel の受理である。五個の有理
first integral、その完全不変体、元の NEW beta quotient との field index、
新しい hypersurface equation はこの receipt の結論ではない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_six_source_review_20261103_v63.py --output work/six-source-independent-review
```

60秒・2.2GB上限。新しい beta 点や高次数 probe は作らない。
