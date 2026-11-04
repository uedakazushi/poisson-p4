# 三つの未決二次 class における m1 の消去

## 結論と範囲

実五生成825関係の Schur 表示について、既知の一次関係37本の span を
定数 column 操作で除いた620本を考える。既知の `variable × linear`
関係および168本の literal d3 の二次関係をすべて零とする三つの quotient
dual を適用すると、残る二次 class は **実 X2 と既知の X1,m0 だけ**で
定まる。canonical cofactor 行列 m の未知一次係数 m1 は不要である。

これは未射影の実 Schur 二次行列の等式ではない。Euler contraction は
係数環に関して線形でなく、その差は既知の d3 二次関係にしか落ちない。
この境界を無視して literal な等式と読み替えてはならない。実 X2 の値を
求めたとは主張せず、全長 T・深さ・nil 三乗の零性も決めない。

## 実 frame と620本の source

基底は `Rbar=QQ[[zA0,...,zA6,rA]]`、最大 ideal は m とする。実 source
frame の分裂を `Rbar^56 -> Rbar^24` の全次数行列 m と、自由32座標に
より固定する。定数 m0 は rank24 である。実 X0 の27正規化独立列を
`Xp` と書く。定数自由32商から得た五成分射影を pi と書くと、

```
m0 Xp0=0,    pi Xp0=0,
ker(m0) intersect ker(pi) = image(Xp0),
```

となる。最後の等号は、自由32と ker(m0) の同型、pi の rank5、および
27列の独立性から従う。pi の元 V56 row ラベルは `[22,23,29,30,55]` で
あり、raw physical 座標への読み替えはしない。

凍結済みの定数27+84 column 操作に続けて、一次 rank37の固定消去を
施す。その620零列の元 source 線形結合を x とする。x0=0 であり、
wedge 定数も零なので m0 x1=0、五生成一次関係も零なので pi x1=0。
従って一意の一次係数ベクトル a1 により

```
x1=Xp0 a1
```

と書ける。a1 は実 X1 と27定数 pivot 行列だけから求まる。builder は
全620本について、この等式を自由32成分だけでなく **全56成分**で確認する。

## m1 と未知 a2 が消える理由

元の全次数 source cycle を保ったまま `z=x-Xp*a1` と置く。a1 は一次
多項式であり、z は m二乗に入る。m は正則形式函数行列なので

```
(m z)2 = m0 C2,       (pi z)2 = pi C2,
C2 = x2 - Xp1*a1.
```

となる。第一の式に m1 が現れない理由は z1=0 であり、任意の m1 を
零に設定したためではない。実27 Schur 消去で使う全係数を
`a=a1+a2+...` とした場合、追加項 `Xp0*a2` は m0 と pi の両方で零に
なる。したがって未知 a2 もこの二次 readback に寄与しない。

## Euler 境界の処理

元の source 列を選んだ contraction を kappa とする。一般には
`kappa(mXp*a1)` は `kappa(mXp)*a1` に等しくない。両者の差の Koszul
boundary は厳密に零で、その差は定数項零で始まる。一次初期形は
`d2_linear` の kernel に属する。

三 sector の一次 Koszul 列は `(u1,...,u7,0)` である。u7 は正則列で、

```
H2 K(u7,0)=H2 K(u7) direct_sum H1 K(u7)=0.
```

従ってこの一次 wedge 差は `d3_linear` の定数係数 boundary である。
同次次数に分ければ、一次係数の cycle を持ち上げるのに定数係数の
G3 だけで足りる。これは g8 の次数が二または三である三 sector 全てに
適用する。CB の `u=z+c*r²` という二次 shift はこの一次複体を変えない。

五生成 Schur 消去でこの差を運ぶと、二次差は
`-Ystar1*U0^(-1)*d3_linear` の span に入る。これらはすべて、既に保存
された168本の実 d3 二次関係に含まれる。そのため三つの quotient dual
Pi3 はこの差を零にする。source の一次関係の選び方の差も、その一次
関係と変数との積として同じ quotient で消える。

## 正確な typed functional

上記の Pi3 は180次元の二次係数空間を三次元へ射影する固定有理行列である。
変数順は `(zA0,...,zA6,rA)`、二次 monomial 順は
`combinations_with_replacement(range(8),2)`、target順は実五成分である。
凍結 rank177 の固定 pivot を再生して Pi3 を取り、全464既知列を零に
することを確認する。新しい rank 検索はしない。

全620 source の実二次 class は、正確に

```
Pi3 [ pi*C2 - Ystar1*U0^(-1)*kappa_linear(m0*C2) ],
C2=x2-Xp1*a1
```

である。ここで Ystar1 と U0 は凍結済みの84単位列の実 frame に従う。
kappa_linear は二次入力に対する一次 Euler readback であり、出力に
作用する Euler を E とすると

```
wi8=(E+1)^(-1) partial_i(v8),
wij=(E+2)^(-1)(partial_i(vj)-partial_j(vi)).
```

整数の正の分母だけを用いる。形式函数に operator を作用させてから
係数を取るのであって、I二乗商の normal form を微分するのではない。

出力は全620 source を元768列へ戻す定数係数、各 a1、実 Xp1、m0、pi、
Ystar1、U0、Pi3 と上記の有限 operator 式を source hash とともに保存する。
**実 X2 callback の値は入力されていない。** 旧 r/g 正則函数の二次微分は、
この functional に依存が残る場合だけ、別の source-order gate になる。

## 純粋な二次 unit-basis 変更の自然性

source56 frame を `P=I+P2+O(m³)` で変え、Q、SyzQ、B と全比較図を
整合的に運ぶ場合、残余 source の左変更は二次で
`-P2*x0=0` となる。m0、pi、Xp0、Xp1 は変わらない。relation frame の
純二次変更による追加項は image(X0) に入るので `Xp0*a2` 型となり、
m0 と pi の両方で消える。contraction や比較 homotopy の選択差は、
前節の既知境界を modulo した class を変えない。

これは一般の **整合的な基底変更の定理**である。native adapter N の
全二次係数が実際にこの変更だけであると判定したわけではない。N2 だけを
SyzQ の式で変え、Q や B を固定する操作にはこの定理を使えない。物理的な
Q2、core pullback の二次変化、または整合的輸送の証明を欠く項は別の実入力
として残る。従って本 packet は native gauge2/source5 の job を自動的に
許可したり、不要と断定したりしない。

## 再現

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_T_quadratic_m1_cancellation_preflight_20261104_v63.py
```

60秒・2GB・32MiB・150,000参照 guard。固定 receipt と既存 X1 のみを
再生し、新しい physical jet、X2 求解、348 scalar 行列、rank探索は行わない。
