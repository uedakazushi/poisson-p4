# 3次・5次の map degree から s を有理回収する

2026-11-04。この小さい独立証明は最終結果式の全係数 identity に依存しない。
受理済み q/B逆式により `r,R,Z0,s` は有理sourceの独立座標である。
`k=QQ(r,R,Z0)` とし、Z1,Z2をk上のs有理函数として見る。

分母はkの単元を除き、それぞれ `J*Fn`、`J*Fn²`。ここで `deg_s J=1`、`deg_s Fn=2`。
分子NA1,NA2のFnでの pseudo-remainderと、NC1,NC2のJでの pseudo-remainderを
独立の sparse Fraction 演算で計算し、四つとも非零であることを全係数で確認する。
FnとJも互いに素である。

Fnの既約性をblack-box factorizationに頼らず確認する。元の
`F=s²+100r(r−1)q²−4rs+4r` はQQ(r)上の非特異conicで、斉次行列式は
`−400r²(r−1)² != 0`。Rの分子・分母の作る直線束の**記号的な中心**は
`q=1/[5(10r−1)]`, `s=10r/(10r−1)`。
そこでのFは `64r²/(10r−1)² != 0`。従ってconicを中心外の点から射影するR写像は
degree2である。qはR,sの有理式なので、sはQQ(r,R)上ちょうどdegree2を持ち、Fnは
その最小二次式となる。これは新しい数値点を用いた検査ではない。

Fn既約、JとFn互いに素、および四つの非零 remainder により、独立なZ0を含む実際の
各分子と分母には共通因子がない。従って

```text
[k(s):k(Z1)]=3,
[k(s):k(Z2)]=5.
```

中間体 `k(Z1,Z2)` に対する `k(s)` の次数は両方の次数を割るため、正確に1。
よって**rを残した実際のgeneric incidence上で**sは
`QQ(r,R,Z0,Z1,Z2)` の有理函数である。
独立PRSのdegree1 remainderが得られ、そのleading coefficientが単元である開集合では
具体的に `s=−C1_constant/C1_linear` と書ける。最終constant resultantの検証は別のgateである。

この結果はrの回収、actual NEW target体のindex1、r射影の次数19、実定義方程式を
証明しない。GB、gcd、resultant、新prime、新点、source26、失敗したPRSの再試行は使わない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_two_shape_coprime_degree_review_20261104_v63.py --output work/v62-review/two-shape-coprime-degrees-independent
```
