# actual Nil(D0) の44成分積・cube operator：独立監査

2026-11-04。受理済み完全複体 `R0^780→R0^44→R0^3` の H1 が
actual `Nil(D0)` であることを固定する。producer はimport・再実行しない。

## 同じ actual source と基底

Qradの全36列が、同じ9つの native `dR9` と4 basis倍であることを独立に追う。
4 basis は `(1,τA,τB,τAτB)`。全3×44 cycle写像について、前36列はQradのodd成分、
後8列はτA行だけがactual β8の負号で、他2行は0であることを確認する。
β8は元C0の7正則関数と元gAの同じ節点、νAも元関数節点20086である。
旧gAを別の近似商に置き換えず、τAとτBを同一視しない。

## 全492節点と普遍積

44＋44個の入力とactual係数を独立記号として扱い、独自の可換多項式演算で
全492節点、44個の双線形出力、44個の三次出力を係数ごと検算する。
source係数関数そのものの新しい展開やIFTはしない。4 basisのbitmask積から
独立に `C=R0[τA,τB]/(τA²−νA,τB²−νB)` の乗法を作る。

cycle `z=(a36,c8)` では `b0=Qrad_even*a`、`B=Σβi ci` とおけば
`Qrad*a=b0+BτA`。保存された積の前36成分は9 blockごとの
`μ(b0+BτA)*a'`、後8成分は `b0*c'+b0'*c` に一致する。
QradのC線形性により、これはactual C0の元の積を与える。

even成分を無条件に `b0b0'+νABB'` と置いていない。
cycle条件を外した一般入力では、その差が正確に
`νA B (Vz')_τA` となり、残3成分も
`[[b0,0,0],[0,b0,νAB],[0,B,b0]] Vz'` と一致する。
したがってcycleの積はcycleである。

## homologyへの降下と正確な三乗gate

完全なH1同型は、cycleを回復したC0代表元のD0-classへ送る。
境界のD0-classは0なので積も0。逆に0-classのcycleは全Xの境界であるから、
積はH1に降りてactual nilidealの乗法そのものとなる。
chain上のa成分が可換・結合的だとは仮定しない。差は完全syzygyにより境界となる。
780列の新しい係数homotopyを求めたとは主張しない。

一入力cubeは `μ((b0+BτA)^2)*a` と
`(3b0²+νAB²)c` に一致し、cycle上でD0の三乗を表す。
標数0の偏極式を独立に確認したので、Nil(D0)^3=0と同値な条件は、
**すべてのexact cycle**のcube44がimXに入ることである。
free44の標準基底をnil生成元やcycleと同一視しない。

結論はsource-bound積作用素と全cycleに関する正確なgate。
nil生成元、三乗の零・非零、H0全体やTor63との同一視、sourcejetや新係数は未主張。

再現：`.venv/bin/python code/ambient_B18_nil_product44_operator_review_20261104_v63.py`
