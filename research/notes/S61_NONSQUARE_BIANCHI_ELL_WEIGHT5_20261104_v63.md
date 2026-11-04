# 全135 Bianchi cofactor の weight5 / ell-divisibility

既存43個の個別変数による除算の否定は、core に依存する線形形式 ell による
除算を扱わない。今回は凍結済み full175+tG6 の真の localweight5 成分20行
だけを、一回60秒・2GBで検査する。

ell は旧 x7 座標の実形式であり、ell_x1=-16CJ/243 は FULLU の既知単元。
係数を normalgraph 上だけで検査すると normalizer cofactor を失うので、
全 linear extended chart を使う：

```
raw45 = T34[:,normal34]*D*M*normal8
      + T34[:,1]*base1 + U34*normalizer34.
```

D は受理済み Tn/T34 adapter、M は v3=W,v6=Z+x1-3x0 を含む可逆定数変換。
normal8 の順は (x0,x1,x2,x3,x4,Z,W,H)。追加 weight1 の base1 と U34 の
34変数を全て独立に保持する。T34 の列2,5は nuisance19,23 の座標単位であり、
それらだけを零にする。Tnormal・Tbase1・U34 の両行が恒等的零であることを
入力から確認する。非線形 normalgraph も source319 も使わない。

全135の各 cofactor はこの43変数で linear または constant。
ell=0 において x1 を既知単元で解き、全 cofactor の全残係数を constraint にする。
F5 の constant cofactor はそのまま残す。F6 と tG6 も135列から除かない。
extra tG6 は真の weight6 なので weight5 成分への寄与は零である。

full constraint layout は5805x20。そのうち非零行の exact QQ rank を、
受理済み同一点 (a,e,d)=(1,2,3) だけで一回計算する。full-column-rank minor が
非零なら、corefield QQ(a,e,d) 上にもこの ell により割れる非零組合せはない。
rank が落ちても generic survivor は主張しない。chart の全weight blockの
同一点 determinant を併せて保存し、generic に可逆な full43 frame を束縛する。
ここで全点の frame-minor 単元性を新しく主張する必要はない。

この否定は上記20行の、全 local 変数について定数・core について有理な係数の
組合せだけを扱う。高 weight、一般 polynomial cofactor、別の線形延長や一般
saturation を尽くさない。新しい uniform kernel 方向・full239には結論しない。
