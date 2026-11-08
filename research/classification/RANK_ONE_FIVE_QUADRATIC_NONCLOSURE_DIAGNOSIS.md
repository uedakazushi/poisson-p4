# M1/M11のquadratic余核非閉鎖：残るのはplane／二line

2026-11-08。元ONE five有限jobの結果を数学的に診断する。旧point／source／caller／reportは不変で、新ambient点、高次jet、degree探索をしない。

全canonical source＋実群作用のFは両点でrank50であり、true normalはM1三本、M11二本。full350 genuine cokernelを使っているため、選んだnormalの任意性又はF方向の取り残しを原因とした偽の余剰ではない。

## M1

保存したnormal順序(t0,t1,t2)で、全quadratic idealは次の三式で生成される。

```
−152*t2*(1377*t0−2394*t1−322*t2)/176769,
−7*t2*(39362922*t0−66490425*t1−8556539*t2)/1785013362,
−4*t2²/27.
```

t2を除いた三つのlinear formの係数行列のdeterminantは−2128/1590921≠0。従って

```
I2=t2*(t0,t1,t2), radical(I2)=(t2).
```

quadratic coneの被約支持は二次元planeである。固定degree四のMacaulay imageは15 monomial中10次元で、t2を持たない五monomialが残ることとも一致する。

## M11

全quadratic idealは一つの式

```
2*(231*t0²−110*t0*t1−70*t1²)/1155
```

で生成される。numeratorのbinary discriminantは76780≠0なので、複素数体上の支持は相異なる二line。固定degree三Macaulay imageは四monomial中二次元である。

## 次に主張出来ること

これらはquadratic initial idealの正dimensionであり、**元のPoisson切片に余剰の正dimension枝が存在する証明ではない**。高い次数の障害がnormal支持を原点へ落とす可能性は未決。逆に接dimension28/27だけから新28/27成分があるとも言えない。

全source25次元像、全priorからの分離、generic軌道及びsource全zero/unit被覆は保持される。今回の有限quadratic路はここで停止する。次は独立なglobal case-bound inventory又は実幾何的親／全支持の構造経路であり、より長いjet、別point又はMacaulay degreeの増加ではない。
