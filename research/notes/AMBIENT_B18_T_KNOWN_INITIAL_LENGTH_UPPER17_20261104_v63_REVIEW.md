# 独立監査：実 cubic 関係による全長上界17

実 source witness71の独立レビューを前提として、新しい canonical cubic
`e4*z5*r²` が真の initial 加群に属することと、40生成の全次数同次計算を
独立に受理する。

独立に構成した GB39 quotient の二つの dual に対して、実 source cubic の
値は `(-251099375/1051776,0)`、canonical monomial の値は `(1,0)` である。
従って前者を非零有理数で正規化して後者を引くと、差は同じ既知 initial
部分加群の次数三部分に入る。単なる数値の非零から未指定式を追加する
操作ではない。source peer は元の768列係数、Xp1*a2の符号、nativeの零-origin
結合範囲、未知 atom二次係数の消去まで確認している。

元39 generatorにこの一項を加え、別の Fraction polynomial engineで最後の
適用可能 divisor を選んで、全140同一target S-pairが零となることを検証
する。作者scriptの import・実行はしない。よってこの40本は既知同次部分
加群の全次数 Gröbner 基底であり、その商の標準項は

- 次数0: 五つの定数 target
- 次数1: e3*r, e4*z5, e4*r
- 次数2: e4*r², e4*z5*r
- 次数3以上: e4*r^d

となる。Hilbert級数は `5+3t+2t²+t³/(1-t)` である。

これは真の全 initial 加群を同定したという主張ではない。その商は gr(T)
へ全射し、同じ実最大 ideal について独立確認済みの m十乗零を使うと

```
10 <= length(T) <= 5+3+2+7 =17.
```

既証の実DVR比較式から mu_D Kは6–13、torsion部の最小生成数は2–9となる。
torsion exponentやH0長の上界ではない。残るpure r³初期classがさらに消える
可能性は排除せず、真のHF3、T全長、全nil三乗の零性は未決のまま保持する。
元の上界24の凍結証明は歴史的な有効結果として変更しない。

再現:

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_T_known_initial_length_upper17_review_20261104_v63.py
```

新しいsource値、jet、rank探索、348行列を使わない小さい固定多項式監査。
