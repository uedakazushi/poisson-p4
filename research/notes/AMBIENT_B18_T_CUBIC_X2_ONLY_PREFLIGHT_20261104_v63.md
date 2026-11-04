# 一次係数 syzygy の三次 class を実 X2 だけで読む

## 結論と限定

実五生成825関係の表示から、固定37本の実一次関係を選ぶ。その一次初期形
に対する119本の一次係数 syzygy が与える **実三次関係 class** は、X2 と
既知 X1,m0 だけから計算できる。X3、m1、m2 の新入力は要らない。

等式は、既知の二次初期関係と変数との積を零とした cubic quotient 内での
等式である。未射影の Schur 三次係数行列の等式でも、全ての三次初期関係を
この119本が尽くすという主張でもない。特に、620本の二次関係の定数係数
kernel に由来する別の三次関係は X3 を要することがあり、本 packet は
それを計算しない。既に受理した原点深さ19と HF(5,3,2) は変更せず、HF3、
全長T、nil三乗の新しい結論は得ていない。

## 全次数 source の設定

実基底は `Rbar=QQ[[zA0,...,zA6,rA]]` である。元768列に対する定数27消去と
定数84消去の後の657本の source cycle を X_j とする。これは全次数の元の
source 列の定数線形結合であり、Schur 係数を任意に零とした列ではない。

```
X_j0=0,       m0 X_j1=0,       pi X_j1=L_j.
```

最初の等式は実 X0 の定数消去、第二は wedge 定数消去、第三は実五成分
Schur 一次係数の式である。builder は選択した37本について、保存された
全56成分 X1 と m0、pi でこれらの一次 interface を確認する。

27本の元の pivot cycle を Xp と書く。既証の実 frame では

```
ker(m0) intersect ker(pi)=image(Xp0),
```

であり、Xp0 の27列は独立である。ここで pi は元の自由32成分から得た
定数五成分射影で、raw physical 座標への読み替えではない。

## 一次 syzygy と二次 pivot 補正

選択37本に対して一次斉次多項式 l_j を選び、`sum l_j L_j=0` とする。
全次数の source cycle

```
x=sum l_j X_j
```

は x0=x1=0 であり、その二次初期形は `x2=sum l_j X_j1` である。
m0 x2=0、pi x2=0 なので、一意の二次係数 a2 により

```
x2=Xp0*a2
```

となる。a2 は実 X1 と定数27 pivot だけから求まる。builder は全119本に
ついて、全56成分でこの等式を確認し、a2 を保存する。これは単なる自由32
成分の kernel 計算に置き換えない。

`z=x-Xp*a2` は全次数の実 source cycle で、m三乗に入る。従って

```
(mz)3=m0*C3,        (pi z)3=pi*C3,
C3=sum l_j X_j2 - Xp1*a2.
```

m1 と m2 が消えるのは z の零次・一次・二次が全て零だからであり、未知
係数を零と代入したからではない。全次数の実27 Schur 補正を a2+a3+... と
書いた場合、未知 `Xp0*a3` も m0 と pi の両方で消える。

## 84単位消去と非線形な Euler contraction

37本の元 source は wedge定数零なので、対応する wedge 関係 W_j は
m に属する。l_j も m に属するため、`sum l_j W_j` は m二乗から始まる。
84 pivot の単位行列を U、下五行を Ystar と書くと、84消去の三次寄与は

```
-(Ystar)_1*(U0)^(-1)*(wedge leading degree2)
```

だけである。U1、Ystar2 が寄与する最小次数は四である。

ただし source の Euler contraction kappa は係数環に関して線形でない。
次の二種類の差を無視してはならない。

```
kappa(sum l_j*mX_j)-sum l_j*kappa(mX_j),
kappa(mXp*a2)-kappa(mXp)*a2.
```

いずれも全次数 Koszul boundary が零で、wedge次数二以上から始まる。
従ってその二次初期形 delta は `d2_linear delta=0` を満たす。各 sector の
一次複体は `(u1,...,u7,0)` の Koszul 複体であり、u1,...,u7 は正則列なので

```
H2 K(u1,...,u7,0)=H2 K(u1,...,u7) + H1 K(u1,...,u7)=0.
```

ゆえに delta は `d3_linear` の一次係数 boundary である。全次数 d3 の
対応する一次倍の leading と一致するため、その五生成三次像は、既知の
実 d3 二次関係と一次変数との積に入る。この span を零とする cubic dual
では両差が消える。CB sector の二次座標 shift `u=z+c*r²` は、入力が三次
から始まるこの二次 wedge readback を変えない。

## 正確な二 channel functional

Pi_cubic は、凍結済み全600三次 monomial の固定 linear rewrite と、既知
178二次関係の変数倍の三次閉包を用いる quotient dual である。五標準
monomial の閉包 rank3 による二 channel を使う。37一次式を Gröbner
基底と仮定せず、保存された全1332 overlap の包含を前提にする。

求める119個の class は、正確に

```
Pi_cubic [ pi*C3 - Ystar1*U0^(-1)*kappa_degree3(m0*C3) ]
```

である。kappa_degree3 の二次出力には正しい output-order Euler 式

```
wi8=(E+1)^(-1) partial_i(v8),
wij=(E+2)^(-1)(partial_i(vj)-partial_j(vi))
```

を使う。元の形式函数へ適用してから leading を取り、348商や I二乗の
normal form を微分しない。実 cycle の u=0 readback が零であることと、
一次 Koszul 完全性により、この leading readback は境界を回復する。

出力は全119 l_j、選択37本の元768列係数、実 Xj1、全119 a2、既知 Xp1 と
上式の型付き coefficient callback を保存する。X2 の値は求めていない。
Pi_cubic の数値配列と source adjoint/pair support は A0 の別 preflight が
担当し、ここでは重複抽出しない。

既存の effective native S2 は各 selected source の定数零結合に対して
有効なので、さらに l_j を掛けた source にも L2 cancellation が成立する。
ただし B11 に入っていない source monomial が必要なら、その代表からは
読み取れない。必要な pair support を先に証明し、未知旧正則商の値や新しい
source jet を勝手に零としない。

## 再現と資源

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_T_cubic_X2_only_preflight_20261104_v63.py
```

全119 full56 pivot 恒等式を既存 first arrays だけで検算する。60秒・2GB・
32MiB・150000参照 guard。新しい X2/X3/source係数、rank探索、GB、348行列は
計算しない。
