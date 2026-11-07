# twisted cubic cone上の全rank零・既約quintic枝：六base三軌道と全linear fiber

2026-11-07、主義務G。
**結論。** twisted cubic cone S⊂P4上で零となり、実退化因子が既約被約quinticである
全rank四Poisson tensorのprojective PGL5飽和像はdimension二十四以下である。
clean、node、normalization/conductorのsmoothnessを仮定しない。

## 1. 全十八源を六baseと十二vertexへ分ける

cone座標をx0,x1,x2,x3,cとし、idealは
`q0=x0x2−x1², q1=x0x3−x1x2, q2=x1x3−x2²`。
全三十rank零源とFrobenius first-jetによるabelian conormal必要条件は
`CUBIC_SCROLL_ABELIAN_CONORMAL_NONDOMINANCE.md` と同じ幾何を使う。
coneの場合、全78×30 linear systemのrankは十二、kernelは十八である。
その完全証明書は `CLASSIFICATION_CUBIC_SCROLL_ABELIAN_CONORMAL_20261107.json` のcases[1]。

全kernelは一意に

`Π=P(x0,...,x3)+∂c∧V(x0,...,x3)`

と分かれる。Pの全sourceは六次元、Vは `Σqi Vi`（Vi任意constant base vector）の十二次元。
元の十八kernel座標におけるbase indexは `[1,6,10,12,13,15]` で、これを(A,B,C,D,E,F)と呼ぶ。
残り十二がvertex sourceである。

全350 projective Jacobiは、c係数からstrict `[P,P]=0` を、
cに依存しない部分から `Euler_base∧[P,V]=0` を要求する。
V²の式は無い。
実際の全350 coefficient行を構成すると、各行はpure base二次又はbase×vertexに分かれ、
混合していないこともexactに確認した。
base=0では全Πがdecomposableでrank二以下なので、rank四解析では全て除ける。

## 2. 全base Poisson支持はlinear三空間とquadric二錐の和

全base二次行のrow rankは九。
保存されたexact rowbasisを非零rational定数倍すると、方程式は正確に

`(2A−F)(3A−F)=0`, `B(3A−F)=0`,

`2AC−2BE+DF=0`, `6AD+2BE−3DF=0`,

`BC=0`, `C(C+D)=0`, `(2A−F)E+2BF=0`,

`12B²−4BE+(C+D)F=0`, `6BD−(C+D)E=0`。              (1)

この全被約支持を次の有限場合分けで解く。
A≠0ならF=2A又は3A。
F=2AではB=0、D=−Cで、Eは任意。
F=3AではE=6B、C=0、AD=4B²となる。
A=0ならF=0、次にBE=0と12B²−4BE=0からB=0。
残りは `C(C+D)=0, (C+D)E=0`。
C+D=0なら前者のlinear支持、そうでなければC=E=0で後者のquadric支持に入る。
逆向きは(1)への代入で成立する。
従って全支持は正確に

`L3={B=0,F=2A,D=−C}`、free(A,C,E)、

`Q2={C=0,E=6B,F=3A,AD=4B²}`

の和である。A=0、nilpotent又はboundaryの例外をgeneric開集合の外へ捨てていない。
nonreduced scheme構造をこの被覆命題の前提にしていない。

## 3. actual SL2作用によりnonzero baseは三orbitだけになる

twisted cubic parametrization `(s³,s²t,st²,t³)` に対するSL2のambient linear fieldsは

`h=3x0∂0+x1∂1−x2∂2−3x3∂3`,

`e=x0∂1+2x1∂2+3x2∂3`, `f=3x1∂0+2x2∂1+x3∂2`。

これらのLie derivativeを全六base basisで直接展開した。
L3とQ2のlinear spanは互いに補空間で、双方が三次元Sym²表現である。
全六×六作用行列とblock分解を
`CLASSIFICATION_CUBIC_CONE_BASE_ORBITS_V2_20261107.json` に保存した。
具体的なbinary quadraticの同一視は

`L3: A*u²+E*u*v−C*v²`、discriminant `E²+4AC`、

`span(Q2): A*u²+4B*u*v+D*v²`、discriminant `16B²−4AD`

である。h-weightはそれぞれ(2,0,−2)で、e,fの二隣接arrowも非零。
これは単なるinvariant polynomialの推測でなく、実際のambient linear actionとの一致である。

複素数体上のnonzero binary quadraticは、相異二根又は平方であり、
SL2とoverall source scalarの作用でそれぞれ一orbitとなる。
よって全nonzero base支持は次の三stratumで尽くされる。

| stratum | dimension | 代表(A,B,C,D,E,F) |
|---|---:|---|
| L3 semisimple | 3 | (0,0,0,0,1,0) |
| L3 nilpotent | 2 | (0,0,1,−1,0,0) |
| Q2 nonzero | 2 | (1,0,0,0,0,3) |

source scalarでbaseをλ倍するとvertex kernel条件もλ倍されるだけなので、
fiber rankはこのSL2×scalar orbit上で一定である。

## 4. 全十二vertex fiberのexact kernel

全350のvertex線形matrixを上記三代表へ制限し、complete kernelとnonzero rank minorをQQで計算した。
その結果は

| stratum | full350×12 rank | vertex kernel | 全fixed affine source上界 |
|---|---:|---:|---:|
| L3 semisimple | 6 | 6 | 3+6=9 |
| L3 nilpotent | 6 | 6 | 2+6=8 |
| Q2 nonzero | 9 | 3 | 2+3=5 |

である。全三百五十行を使用し、選んだsubmatrixのrankだけからnullityを推測していない。
各complete kernelの代入零とrank minorを共に保存したので、rankの上下界が一致する。
各代表が全orbitを代表する証明は§3であり、三つの任意点による標本検証ではない。

## 5. 全PGL像は二十四以下

coneを保存する群はdimension八である。
binary cubicのbase PGL2三次元、vertex相対scaling一次元、
`c→a*c+Σb_i x_i` の四translationが忠実に作用する。
この群は全rank零・abelian-conormal sourceを保存する。
全nonzero base源のdimensionは§4で九以下なので、projective PGL5像は

`24+9−1−8=24`

以下となる。各nilpotent stratumではさらに二十三、二十以下である。
base=0源はrank二以下であり、rank四点は別に漏れない。
これでtwisted cubic cone上の全rank零・既約quintic枝のG非優越が従う。

## 6. 有限計算・旧負結果・初回ラベルの訂正

算術producerは `code/classification_cubic_cone_base_geometry.py`。
base九式と全linear vertex matrixは `CLASSIFICATION_CUBIC_CONE_BASE_GEOMETRY_20261107.json`、
正しい三orbit証明書は上記**ORBITS_V2**である。
V2の全作用・三kernel算術は約一秒で完了し、GB又はhigher-power計算を使わない。

初回ORBITS reportは、exact invariantがE²+4ACであるにもかかわらずA=1をsemisimpleと誤ラベルした。
そのreportはnilpotent二点を調べたarithmetical diagnosticとして不変保存し、
全stratum被覆のclaimを無効とするERRATUM reportを添えた。
V2はE=1をsemisimpleとしてinvariant非零をassertした後に全kernelを再計算した。
旧reportのuniversal booleanは証明に使わない。

旧三十sourceのrowspace負結果、GB timeout、平方test、及び新十八source局所化timeoutは全て保存している。
今回の閉鎖は、新しいFrobenius線形必要条件と全base support/actual orbit/linear fiberに依存し、
それらの負結果を成功に読み替えたものではない。
