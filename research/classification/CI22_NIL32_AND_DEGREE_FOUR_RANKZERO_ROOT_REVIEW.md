# degree四・rank零の既約quintic枝は全空：統合側の独立確認

2026-11-08。PRIMARY `E_reduced_nonnormal`。証明読解と手計算による確認であり、producerの再実行を独立算術とは数えない。

## 六源の全点での必要因子

`CI22_NIL32_RANKZERO_ACTUAL_FACTOR.md` の座標・符号を採用する。
全源は `P04=f, P13=-f/2, P03=-g, P12=g, P01=-k`、
`f=a0*q1+a1*q2, g=b0*q1+b1*q2, k=c0*q1+c1*q2` である。
`q1=2x0*x4-2x1*x3+x2², q2=2x2*x4-x3²`。
全六係数を保ち、Jacobiの四添字1234を直接展開すると

```
(E wedge J)1234
 = -x4*((a0*x2+b0*x3)*f+(a0*x3+2*b0*x4)*g).
```

この式の `x0*x2*x4²` と `x0*x4³` の係数は、それぞれ `-2*a0²` と `-4*b0²`。
したがって全複素Poisson閉点で a0=b0=0 である。被約台だけの推論であり、元idealに a0,b0 自体が属するとはしない。

同じ五つの括弧のEuler Pfaffianを独立に外積すると

```
F = -x2*f²/2-x3*f*g-x4*g²
  = -q2²*(a1²*x2/2+a1*b1*x3+b1²*x4).
```

`k` の項は対となる四添字Pfaffianが零なので全て消える。
非零rank四点は正確に `(a1,b1) != (0,0)` の側にあり、rank三quadric q2 の接hyperplaneと二重q2という実因子型を持つ。
全H+2rank3Q支持定理はこの配置も全点として含む。
残る `(a1,b1)=(0,0)` は decomposable tensor で一般rank二以下。
従って全rank四点は既知閉和へ入り、既約被約quinticという元の枝は全空である。
未解の残Jacobiを自由と仮定せず、必要二係数と全実Pfだけで足りる。

## 全degree四支持への接続

`DEGREE_FOUR_RANK_ZERO_SURFACE_EXHAUSTION.md` の全centre還元と別作者reviewを適用する。
integral・非退化degree四曲面は、CI22、projected Veronese、S(1,3)、S(2,2)、S(0,4)の全外点射影で尽くされる。
nonCI枝は既存の全源証明で空。全member退化CI22も既存のKronecker全型証明で空。
非退化memberを持つCI22はself-adjoint Jordan分類によって尽くされる。
二固有値以上はcommon-imageによる実一次因子、単固有値5と4+1は実二次因子で空であり、
残っていた3+2だけを上の二Jacobi係数が閉じる。
よって **実Dが既約被約quintic、Sがintegral非退化degree四、πのS上restrictionが零という全枝は空**。
古い像24以下というG結論を全点へ読み替えたのではなく、その唯一の残源に新しい必要因子定理を接続した。

rank二restriction、degree五・六、plane支持、及び全Eの残源はこの結論に含めない。
sourceの全Euler/Koszul還元と全標準形は既受理証明に依存し、今回改めて全分類を独力で算術再生成したとはしない。
新exact全350×21/126×21とinline実装はproducer report、入力bytesと旧V2は13入力manifestに保存されている。
