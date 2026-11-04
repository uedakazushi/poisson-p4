# Laurent degree width による共通 content の追加 screen

2026-11-04。29段後の凍結残余を読み取るだけで、多項式係数の積・差・除算、
gcd、pivot、級数、新しい source 計算を行わない。

各非零 entry P に対し、変数ごとの Laurent exponent の幅
width_v(P)=max exponent_v−min exponent_v を数える。共通 Laurent 因子 g を
各変数の最小冪で正規化すれば、g と各 complementary factor は各変数の最小冪0の
多項式になる。整域上の各変数次数は積で加法的だから

```
degree_v(g) ≤ min_(slice entries P) width_v(P)
```

が成り立つ。五つの上界が全て0なら、共通非単項式因子は存在しない。
一つの単項式 entry がある場合だけでなく、変数ごとに別々のentryが幅0を持つ
場合にも使える。全最小値と、それを実現する元source成分番号を保存する。

正の上界が残ったsliceについては、正規化gのdense支持数の上界
product_v(degree_bound_v+1)を計上する。これはgの存在やその支持を計算した
結果ではない。gcdやprimitive quotientの中間係数上界は別に必要であり、
新しいgcd workerは実行していない。

行列、全五parameter、全制約と既存unsafe forecastは不変である。
