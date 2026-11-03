# Actual135からの小さいlazy quartic coefficient oracle

2026-11-03。全四次122×64判定の33genuine quarticへ直接アクセスするための実係数DAGを構成する。大きなparameter polynomial matrixは展開しない。既存CloudのS61 sourceとweighted chartだけを再利用し、GB、新しい実core点・有限体、追加の未知unitを使わない。四次rank64やambient除去の証明ではない。

## 1. 元135式と既知unit

元135式を現在のbetaとfixed q-coordinateでspecializeし、local45 coordinatesにおいて

`F(x)=J*x+B(x,x)`

と分ける。Bは元のquadratic polynomialをpolarizeした対称bilinear mapで、異なる二変数の係数を1/2ずつ両方向へ置く。元sourceの行134（G6）も含む全135行のJとBを保持する。

既存weighted chartの8normal columnsをT、選択34行のright inverseをUとする。厳密なcore arithmeticにより、全135行についてJT=0、選択34行についてJ_selected*U=I34を再検証する。nuisance local coordinates19,23はTとUで恒等的に零。T,UとJ/Bから使う全rational atomsの分母は、非零QQ定数と既知a,F,B,C因子だけであることを検査する。実際のsource条件はaFBC!=0で、DeltaHやDeltaQをこのoracle生成には使わない。

## 2. 正の重みによるexact finite graph

normal v0,...,v7の重みは `[1,1,1,1,1,2,1,1]`、H=v5。E=T*vとし、選択equation coordinateをwとする。x=E+U*wと書くと選択式は

`w+B_selected(E+U*w,E+U*w)=0`。

ordinary polynomial degree別に

`w2=-B_selected(E,E)`,
`w3=-2*B_selected(E,U*w2)`,
`w4=-B_selected(U*w2,U*w2)-2*B_selected(E,U*w3)`

となる。選択34行の最大weightは4であり、全normal variableのweightは正。したがってordinary degreeが5以上のmonomialはこのweight範囲に存在しない。各local variableとselected equation weightが厳密に保たれることも検査する。

従ってw=w2+w3+w4は、単なる4jetではなく、全選択34式に対する**この正weight chartのexact polynomial graph**である。重みを上げる帰納法でも、weight1のwは0、各次のwは既に求めた低weightだけから一意に決まる。ここでformal jetの一致から全環の一致を推測していない。

## 3. Residual rootsと33 quartics

P=J*Uと書き、raw bilinear expressionsのordinary degree2/3/4成分を

`b2=B(E,E)`, `b3=2*B(E,U*w2)`,
`b4=B(U*w2,U*w2)+2*B(E,U*w3)`

とする。各成分のresidualは `b_k-P*(b_k)_selected`。これをDAGで表し、選択34行のrootはcanonical cancellationで零となることを検査する。

全135行の元J/Bを保存し、weight<=4のresidual rootsを指定する。選択34行以外のweight4行は正確に33本で、全てのoriginal row indexとrootが入る。weight>4行を零と主張しない。元135のwhole idealを四次以下と同一視することもない。

DAGはvar、core constant、linear combination、productだけを持つ。各nodeのweighted degreeとordinary degreeを追跡し、weight>4のproductを除く。core atomsはQQ(a,e,d)の既知unit分母を持つ実係数であり、parameter polynomialの巨大展開を行わない。

## 4. Coefficient API

```python
from s61_nonsquare_lazy_quartic_oracle_20261103_v63 import Oracle
oracle, data = Oracle.load('reports/S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63/oracle.json')
row = data['genuine_quartic_original_rows'].index(77)
root = data['genuine_quartic_roots'][row]
node = oracle.coefficient(root, (0,0,0,0,0,2,0,0))
```

返却nodeは同じoracle上のweight0 arithmetic circuitであり、この例ではH²係数。`coefficient`はlight monomialの有限convolutionとmemoizationを使い、coreを展開しない。必要な場合だけ `oracle.core_value(node)` を明示的に呼ぶとQQ(a,e,d)で展開される。これはoptionalであり、全係数を一度に展開するAPIではない。

source normal variablesの順は従来のv0,...,v7と同じ。122 circuitで使うlight順は `[v0,v1,v2,v3,v4,v6,v7]`、H=v5である。actual quarticを正規化genuine cubic basisへ置き換えず、元row番号のままqueryできる。

## 5. 検証範囲と次の利用

producerはsource/core identities、全weighted compatibility、選択34 rootsの厳密な零、全33quartic roots、serialized DAG roundtripを検証する。coefficient抽出は別の小さいuniversal polynomial fixtureの全weight4 monomialに対して直接展開と照合する。このfixtureは実core点の検査ではない。

次はこのoracleから、既存122 circuitの必要な少数係数や因子だけをqueryする。heavy7/pure146 chartを通す実122×64のrank64は未証明であり、coefficient oracleの完成と混同しない。

再現は `python code/s61_nonsquare_lazy_quartic_oracle_20261103_v63.py --output <new-output>`。DAG packetに8MiB上限、node数に100000上限を設ける。新しいツール導入や設定変更は行わない。
