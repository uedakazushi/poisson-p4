# reduced R・非共点型m三の九つの全higher源：最弱線形上界と因子

2026-11-08、PRIMARY G、root producer。独立幾何・算術review前。
対象は被約既約な実五次退化因子、rank二restriction plane、normal位数三である。
零三次曲線Rは接線conic+line又は非共点triangleとする。
横断conic+lineは別幾何証明で空。concurrent又は非被約R、全点E帰属は含めない。

## 全源への入口

`PLANE_CONIC_LINE_M3_ROOT_COVERS_AND_EXCLUSION.md` はconic+lineのsplit、C2、C3と
repeated全root型を尽くし、横断型を空、接線型を四つの全first-symbol軌道S/N/D/Mへ還元する。
`PLANE_TRIANGLE_M3_ROOT_COVERS_AND_FIRST_SYMBOLS.md` はtriangleの全root coversを五つのtoric軌道へ還元する。
これらの幾何的尽くしは計算とは別の依存であり、独立reviewを要する。
literal全P0とleading targetはそれぞれtwo first-symbol JSON inputsに固定する。
leadingの係数を一点に固定せず、全parallel vector spaceへの所属だけを課す。

## 一回の全九源計算

`jobs/plane_reduced_m3_nine_fullhigher.json` のgate PASS後、既general runnerにisolated modeを加え、
有限CPU230秒/RAM2048MiB内で一回だけ実行した。
旧六FunctionDefと非dispatch module statementsは不変、旧df9c全文及び当時全36契約は
`inputs/plane_reduced_m3_nine_frozen_contracts_20261108.json` に不変保存してから追記した。
既存のsmooth cubic計算をやり直すものではない。

全raw一次normal係数36/Euler六を忠実三十代表にし、全二次normal係数九を残す。
元projective Schoutenの全350行からJ1（350×30）及びJ2（350×9）を直接再構成した。
raw strict Jacobi零は仮定しない。特にMも元projective全式を用いる。
actual Euler-Pfaffianのnormal次数三はP0/P1 crossに線形であり、P1²とP0/P2は次数四以上なので、
全24×30行と全parallel target annihilatorをJ1へ追加して完全核K1を得る。
固定P1上のP2解は空又はker J2のtranslateの部分集合である。
高いJacobi条件を省略することは必要supersourceを大きくするだけで、上界に有効である。

固定planeで全first-symbol軌道の次元Bも、全150×39行から確認した。
三十九列はEuler二十五、base/normal Levi十三、tensor scale一である。
Euler rank二十四を引いた値がBで、全normalizerとliteral軌道への幾何的還元がその適用を保証する。
plane選択六とtensor射影化一を戻すと、全像の上界はB+dimK1+dimK2+5となる。

| 全profile | B | dim K1 | dim K2 | 全像上界 | さらに必要な実因子 |
|---|---:|---:|---:|---:|---|
| S | 7 | 8 | 3 | 23 | u |
| N | 9 | 8 | 3 | 25 | u |
| D | 9 | 9 | 3 | 26 | u |
| M | 10 | 7 | 0 | 22 | 因子判定を用いない |
| SF_C3_UNUSED | 9 | 8 | 1 | 23 | 因子判定を用いない |
| TRIPLE_HX2 | 9 | 8 | 3 | 25 | u |
| TRIPLE_HY2 | 9 | 6 | 1 | 21 | u |
| TRIPLE_HYZ | 9 | 9 | 3 | 26 | u |
| TRIPLE_HXY | 9 | 7 | 1 | 22 | u |

各profileの全P0はu=0を保存する。P1のu-rowのu非整除部分はv²Cだけで、
Euler六gaugeでは変わらず、P2にはnormal columnが無い。
七profileでは完全K1の全三bad-C座標が零だから、元全tensorのu-rowがu倍となる。
u=0上は三次元射影接空間内のbivectorとなりrank高々二、従って実Dがuで割れる。
実Dは非零既約五次という仮定と矛盾し、これら七つの全源は空である。
単一kernel vector又は特定P1上だけの因子ではない。
残MとSF_C3_UNUSEDは全像それぞれ22、23以下で、未知25–27優越成分の一般元にはなれない。

## 証明書と限界

`CLASSIFICATION_PLANE_REDUCED_M3_NINE_FULLHIGHER_20261108.json` に全原行、完全核、全行消滅、
非零rank minor、target/annihilator、全Euler/Levi/scale行、実行producer全文と実測資源を保存した。
全九profileの判定は成功した。quadratic compatibility、GB、jet又は追加profileの計算はしていない。
幾何的全coverと上の算術の独立受理が揃った時点で、全非concurrent reducedR/m三のG非優越を登録する。
二つの非空supersourceの全点E帰属・実解の存在・ambient成分性はこの上界から従わない。
