# Actual両zero-odd idealの和はlength3：T3は有限長

2026-11-04。新しいjet計算を行わない、combined105証明書の正確な係数からの帰結。

R0を共通physical九変数完備局所環とする。actual CAの七つのzero odd座標をz7、
残り二つのactual座標をr=j_CA46,s=j_CA95とする。
CA zero odd idealは`I0_CA=(z7,g_CA(r,s))`で、g_CAのorderは3である。
同じphysical baseのCB zero odd8本をz=0に制限しB_i(r,s)と書く。
初回combined105系とその三次証明書が、八つ全てに定数・一次項が無いことを与えた。

三つのactual restriction B18,B60,B67について、(r²,rs,s²)係数行列は

```
[-95/72,    -115/36,  -5/3]
[-1385/288, -875/72,  -145/24]
[-65/72,    -10/9,     0]
```

であり、その正確なdeterminantは **625/768**。
従ってF=(全八B_i)はm=(r,s)に対して

```
F ⊂ m²,      F+m³=m².
```

有限生成module M=m²/FにNakayamaを適用する。
上式はM=mMを意味するのでM=0、従って **F=m²が全次数で正確に成立する**。
これは有限jetから一般のmembershipを推測する手続きではなく、有限生成idealに対する
単位二次係数行列とNakayamaによる証明である。
g_CA∈m³⊂m²なので、actual idealの和の商は

```
R0/(I0_CA+I0_CB) ≅ QQ[[r,s]]/(r,s)².
```

従ってlength3、Hilbert function(1,2)である。
二つの個別normal formのabstract座標を同一視したのではなく、literal105系が同じphysical xを
共有することが本帰結の要点である。

先行のsymmetric-annihilator証明からI0_CAとI0_CBは両方T3を消す。
従ってT3はこのlength3環上の、十二generator以下の有限moduleであり、

```
length(T3) <= 36,
m_R0²*T3=0.
```

この長さ上界は前packetのparameter閉fiber上界とは異なり、T3自体の有限長上界である。
ただしT3が非零とはまだ分からない。非零を別のactual membership obstructionで示せば、
T3は全て有限台torsionなので全横断Dのdepth0、cone depth19が従う。
今回そのnonmembership、実長、socle次元、全Dのnilは計算していない。
T3が零だった場合も、他のT0,T1,T2の判定は依然必要である。
ambient二行matrixをI0で減らすTorを無視した操作は行わない。

専用source `code/ambient_B18_T3_double_annihilator_length3_20261104_v63.py`は
凍結済みjoint-degree3-certificateを読み、三行列のdeterminantと全lower-coefficient零を確認する。
新たなjet、点、GB、source再生成はない。
