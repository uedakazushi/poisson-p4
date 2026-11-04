# B18：二次正規化による nil 指数上界12の独立証明監査

2026-11-04。実交叉 `D=CA∩CB` の既定 B18 点の完成環について、
冪零指数の範囲 **3以上12以下** を受理する。深さ19の結果、Tの全長、
比較核の generic rank や支持分類を上界の証明に使わない。

## 任意の基底変換での核 ideal

包含 `A⊂B=A+Aτ` と `τ²=sτ+ν`、`s,ν∈A` を仮定する。
任意の A 代数 M は、有限とは限らない変数集合に対して `A[X]/J` と書ける。
核の代表 k が `JB[X]` に入るなら `k=j0+j1τ`、`j0,j1∈J` と書ける。
`a=k−j0` とおくと、B[X] 内で

`a²−s j1 a−ν j1²=0`。

A[X]→B[X] は包含なので同じ式は A[X] でも成立する。従って `k²∈J`。
2が可逆なら偏極式によって核の任意の2元の積も0になる。
個々の元の平方零だけで止めず、核 ideal 全体の平方零を得ている。

独立 checker は6変数の自由可換多項式環で、次の恒等式と偏極を係数ごと確認する。

`a²−s j1 a−ν j1² = (a−j1τ)(a+j1τ−s j1)+j1²(τ²−sτ−ν)`。

B の A 平坦性、基底変換後の単射性、A上の基底 `(1,τ)` の独立性、
ν・判別式・非単元 source パラメータの逆数はいずれも使わない。

## 同じ actual source の2段階への適用

受理済み同一 physical R0 上の source は
`Ai=R0+Iiτi⊂Bi=R0+R0τi`、`τi²=νi` である。τA と τB を同一視しない。
R0 は完備 Noether 環、Ii は有限生成なので Ai、Bi は R0 上有限。
共通 zero ambient 上の tensor は R0 上の有限 tensor の商であり、既に完備である。
従って actual completed tensor の矢印は、単なる形式的類比ではなく
それぞれ普通の有限基底変換の矢印である。

`D0 → D1=BA tensor AB → E4=BA tensor BB`

各段階の核は平方零。合成核 K の第1段階への像の平方が0だから
`K²⊂ker(D0→D1)`、従って `K⁴=0`。Torを落としたり、中間像を環全体と
同一視したりする必要はない。

既受理の whole E4（coreだけでなく embedded localization を含む）の
nil ideal は三乗零である。環準同型は nil 元を nil 元へ送るので
`Nil(D0)³⊂K`、従って `Nil(D0)^12=0`。

## 全完成環と下界

独立検証済みの同じ weighted model では全正次数の係数写像が stable tail に
単射であり、その tail の nil ideal は平方零。よって nil homogeneous 2元の
積は、総重みが正なら0になる。零次数と正次数との混合積も含む。
閉じた torus-stable nil ideal の有限 Vandermonde 射影と完備性により、
この結論は任意の形式級数にも成立し、q≥2では nil ideal のq乗は零次数のq乗と等しい。
有限 w-polynomial model を先に局所化するという受理済み Laurent scope を守る。
実 G19 power-series 因子の復元と literal scalar の除去は冪零指数を変えない。

Qb の実 family は P0 を通り、その一般点の whole D で nil 二乗が非零。
coherent nil-square ideal の閉支持と faithful completion により P0 でも非零なので、
下界は3である。従って全 actual B18 germ の指数は `[3,12]` に入る。

三乗の零・非零、上界8、正確な指数、全大域交叉や全 Poisson scheme への
主張はしない。深さ19から nil 三乗非零を推論することもない。

再現：`.venv/bin/python code/ambient_B18_quadratic_normalizer_nil_upper12_review_20261104_v63.py`
