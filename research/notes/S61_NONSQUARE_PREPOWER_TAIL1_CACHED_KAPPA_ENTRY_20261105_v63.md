# 追加tail1の一回計算：T1は7595項の非零多項式

2026-11-05。親の一回GO、独立ROOT_SOURCE_REVIEWとROOT_CONTROL_REVIEW後、
immutable guard_v2 workerを一回だけ実行した。他の七列は評価していない。

実式はcurrent15四行`[2,4,5,8]`・original32列1の

```
T1=B*Σκ_i*M[row_i,column1].
```

149-node sourceとliteral9268係数、保存済みκ四個4920係数だけを使う。
旧checkpoint/Q2/572sourceをscalar workerへ読み込んでいない。
H2/H3・anchor22・B clearing・全tailのframeは先行source familyに固定済みである。

4.547秒で三passを完了し、**T1は7595項で非零**となった。
peak owned58986、整数積24571908対、入力を含むRAW上界2269813 bytes。
source359-bit／91676 cells／9517864 bytes／45秒soft・60秒hard／2GBを維持した。
forward後に逆順で全四積を引きliteral零を確認し、forward再現後に一つの出力を保存した。
四source列と最終T1を全件readbackした。部分prefixは発生していない。

出力は`integer-column4-and-result.txt`のJ1KV001で、各多項式は一つの正整数Dと
packed key、整数係数cの全列により`Σ(c/D)*a^i e^j d^k`を厳密に表す。
指数は`key=i*1048576+j*1024+k`。最後のnode155がT1である。
QQ分数係数を重複JSON展開しない。保存した三pass検算はself-consistencyであり、
実source式の独立評価と同一視しない。

既存のsource恒等式は`A*D15_1=-T1*D12`。従って実式の独立監査後には
FULLUの`D(T1*D12)`も元full cubic57の十分開に加えられる。
T=0で新しく点を覆う、TとT1が互いに素、実bad Fitting idealと等しいとはまだ言わない。
T=0/D12=0の補集合、既存Theta/P80/P100/Bbad13の別被覆、全FULLU57・全quartic239・
一般五点尽くしは区別して保持する。次は保存済み係数のsource-only支持比較と、
補集合改善を証明するための有限計算gateを検討する。
