# T3の両側zero-odd annihilatorと、全186 DAGの原点override

2026-11-04。実際のell-colonはまだ未計算。既凍結DAGを変更せず、独立に決定した四つの
g商原点値0を適用し、次の実データと新しいannihilator還元を得た。

## 原点値の別記録

`AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63`の四つの値をnode
330164,330289,330524,330594に適用した。coordinate-division nodeには既存の証明済み
原点値を使い、全339669 scalar nodeを有理数で評価した。
全unit inverseの分母が非零であることも再確認した。
Q[3]の2034列とtop12列は、原点で全て(0,0)である。
これはambient二行matrixの原点評価であり、T3のfiberが零であることを意味しない。
既知のTor問題を回避する証明には使っていない。

## 新しい厳密なannihilator

actual Dは同じphysical base R0[[w]]上の正次数付き環であり、wは共通のphysical row140。
CAをambientとした既存presentationでは、degree3以降の共通module Hからの関係が
degree4で安定する。従って

```
T3 = ker(w:D3 -> D4)
```

はambientの選択に依存しない。
CBをambientとしても、actual CAの完全kernel99本はweight≤2、残りのphysical graph87本は
weight≤3である。係数を持ち上げる関数はweight0で、これは既受理の同じminimal20座標と
zero/positive source分解から従う。CAとCBの二つのtrace basisを同一視する必要はない。
従って同じdegree4安定化がCB側にも成立する。

CB側でH_B/C0_Bはpositive odd coefficient二本を掛けたtrace odd classで生成される。
その各classをzero odd f_Biで掛けると

```
f_Bi * ((j_Be*tau_B)*c) = j_Be*((f_Bi*tau_B)*c)
```

となり右辺はdegree3のC0_B関係である。
従って、既受理のI0_CA*T3=0に加えて **I0_CB*T3=0** も成立する。
これは両traceを混ぜた二行matrixの計算や、個別normal formの無断座標同一視ではない。
したがってT3は新たにactual環

```
R0/(I0_CA+I0_CB)
```

の有限moduleとして扱える。

## 実際の一次データ：この段階では有限台にならない

既保存のactual CA/CB各10×9 Jacobianから、zero odd各8本の計16行を同じphysical R0で比較した。
その合同行列のrankは **7** のままである。
actual CAのF9=(z7,r,s)への逆Jacobianを掛けると、CBの八つのzero odd functionの
z=0平面上の一次部分は全て零となる。
従って新しいannihilatorが一次でNoether parameterを与えるという期待は失敗した。
二つのidealの一致、共通curveの存在、有限台性のいずれも一次情報だけからは分からない。

次の具体的入力は、一つのCAcurve57系へ同じphysical xを通じCBsource48系を接続した
105-unit系と、その上の八つのCB zero odd restrictionである。
これらとactual g_CAの初期成分を比較すれば、共通接線や有限台性への新しい実証明を試せる。
今回のpacketはその高次restrictionを計算したものではない。

## 証明書と範囲

`code/ambient_B18_T3_double_annihilator_preflight_20261104_v63.py`、
`reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/`に、全原点値、
全unit分母、全Q3/top residue、両actual zero idealの一次行列、CA逆Jacobianと制限を保存。
一回0.666498秒、最大RSS151676KiB。GB、新しいjet、点、source再生成、p逆元は無い。
T3の零性、socle、ell-colon、全Dのdepth/nil構造は未決定である。
