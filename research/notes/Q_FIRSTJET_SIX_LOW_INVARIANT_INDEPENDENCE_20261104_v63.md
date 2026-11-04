# F4,p1,p2,p3,p4,N12 の六つの代数的独立性

2026-11-04。既凍結の**同じ beta 点 index0**のみを使用する。新しい点、native
probe、高次数の不変式評価は無い。受理済みの `(Inv8 五個,m)` の6×16微分と
非零6 minor、Inv8→pencil5 の全係数変換を入力として、一つの exact QQ minor
を保存する。

重要な恒等式は全 W16 上で

```
p0=detC=(tr(D²)/2)²−4detD=Delta8,
F4²=400p0+1600p4,
N12=2Delta8*m=2p0*m.
```

最初は受理済み universal f3 と4×4 Newton identity、二番目は受理済みの
全 pencil polynomial identity、三番目は m の定義である。点で値が等しい
だけの式を微分しているのではない。

同じ点では p0 と F4 は共に非零である。従って上の恒等式を微分し

```
dF4=(200/F4)(dp0+4dp4),
dN12=2m*dp0+2p0*dm
```

を使う。`(p0,p1,p2,p3,p4,m)` から `(F4,p1,p2,p3,p4,N12)` への微分変換の
行列式は `400p0/F4`、その前の Inv8→pencil5 も可逆。元の同じ6列の minor
にこれらを掛けた非零 QQ determinant を全36成分・全6×16微分と共に保存する。

従って**この六つの特定の**普遍不変式は代数的独立である。その六変数で
非零な polynomial は、全 beta 上でも非零な polynomial になる。quotient の
次元6や五個の八次不変式の線形独立性だけから推定した結論ではない。

この packet 自体は NEW 上の degree24 relation を計算せず、その存在も
主張しない。後続で非零 polynomial が source 上で厳密に零と証明された際の
「普遍恒等式ではない」という gate のみを供給する。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_firstjet_six_low_invariant_independence_20261104_v63.py --output work/six-low-invariant-exact-minor
```
