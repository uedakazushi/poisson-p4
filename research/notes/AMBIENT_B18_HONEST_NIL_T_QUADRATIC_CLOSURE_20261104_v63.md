# 実 T_nil の既知二次閉包と完全814 source projection

2026-11-04。実35本の一次関係に八変数を掛けた280本の既知二次初期関係は
**rank174**で、180次元空間のquotientは**六次元**と得た。
一回の固定QQ計算は3.30秒・約20MB、33,271 QQ slots・248KBで完了した。
新しいX3、m2、source/gauge jetは一切計算しない。既知閉包のrankを真のHF2と同一視しない。

完全minimal表示M5×849は原点零で一次rank35。任意形式係数a=a0+a1+...について、
二次初期関係となる条件はM1 a0=0で、その二次値はM2 a0+M1 a1である。
従って全二次初期関係は、上の280本と、M1の完全constant kernel814本に対する実M2値に
尽くされる。全814 source組合せ、既知閉包の完全rank証明、残るquotient dualを保存する。
このsource theoremは未計算の二次関係を見落とさないためのものである。

## d3の二次寄与はまだ追加しない

全168 d3の一次値は既に35本の一次spanへ含まれる。二次initialへ追加するには
一次を消すsource組合せが必要で、literal d3二次項を単独で関係と称してはいけない。
さらに81 wedge単位対の消去後には、残る五行C0=0でも

  d3_min,2=d3_keep,2−C1 Wunit0^-1 d3_pivot,1

という補正が残る。C1は81pivot X列の一次値であり、既に求めた最終681列の一次値とは違う。
最終三wedgeは完全線形なsector0なのでliteral d3_keep,2は零だが、補正は未知である。
この未知補正を零に設定しない。

残るquotient dualを全814 sourceの実M2へ適用すれば、必要なscalar adjointを厳密に指定できる。
一般にはh2がsector0 row7の(mX)3を読み、X0c=0ならm0 X3c+m1 X2c+m2 X1cが残る。
保存済みm1は保持し、B11に無いpure-r三次由来m2やpivot補正が消えるかは別のsource証明を要する。
全M1零sourceとhomotopy境界を使う低次数化の可能性はあるが、本packetでは断定しない。
完全GBを計算していないので、有限長上界やannihilatorも新たに主張しない。

残る六座標のラベルはtarget0のz5二乗・z5r・r二乗、target1のr二乗、
target3のr二乗、target4のr二乗である。実quotient dualは単なる座標選択ではなく、
保存された疎な六行を用いる。これらを実M2へ適用する前に、source感度とcache範囲を検査する。
