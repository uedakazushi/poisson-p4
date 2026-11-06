# 冪零(3,1)型・長い鎖の定数項：全六源と二つの22次元成分

2026-11-07。主課題 E_multiple_hyperplane。
modular 場を `X=x2∂x1+x3∂x2+∂x3` とする全 affine 二次 Poisson 構造を扱う。
入口は `CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json` case1 の
全11源と全 Jacobi である。全90線形係数還元は当該証明書を入力とし、ここで再実行しない。
本ノートはその全残条件を解き、二源を既存成分へ認識し、さらに二つの新全成分を得る。
r=−1 の二源は未被覆のまま明記する。

## 1. 全六源の導出

入口の11変数を `n0,...,n10`、r=n10 とする。
n2,n3,n7,n9 は全ての枝で自由であり、全 Jacobi の残条件は

\[
\begin{gathered}
 n_5(r+1)=0,\quad n_1(3r+1)=0,\quad
 2n_1n_9+(2r+1)n_6=0,\\
 (r+1)n_0-n_3n_8+n_5n_7+n_6n_9=0,\quad
 n_8(3r-2)=0,\quad (2r-3)n_4-5n_5n_9=0.                 \tag{1}
\end{gathered}

r が `−1,−1/3,−1/2,2/3,3/2` の外なら、(1)は順に
n5=n1=n8=n6=n4=n0=0 を強制し、実 Pfaffian が零となる。
従って rank四ではこの五値だけを残せばよい。
相異なる五値を代入すると、除算で特別な枝を失うことなく次の全六源を得る。

| 名称 | r | 五つの自由変数 | 他の変数への全代入 | rank四条件 |
|---|---|---|---|---|
| A | −1 | n0,n2,n3,n7,n9 | n1=n4=n5=n6=n8=0 | n0≠0 |
| B | −1 | n0,n2,n3,n5,n9 | n1=n6=n7=n8=0,n4=−n5n9 | (n0,n5)≠(0,0) |
| C | −1/3 | n1,n2,n3,n7,n9 | n4=n5=n8=0,n6=−6n1n9,n0=9n1n9² | n1≠0 |
| D | −1/2 | n2,n3,n6,n7,n9 | n1=n4=n5=n8=0,n0=−2n6n9 | n6≠0 |
| E | 2/3 | n2,n3,n7,n8,n9 | n1=n4=n5=n6=0,n0=3n3n8/5 | n8≠0 |
| F | 3/2 | n2,n3,n4,n7,n9 | n0=n1=n5=n6=n8=0 | n4≠0 |

r=−1 では残る積条件が n5*n7=0 であり、A と B はその両枝である。
交わりを除外せず、片方を他方の稠密部分とみなしていない。
各源の指定された rank四開は既約である。

座標を `(x,y,z,w)=(x1,x2,x3,x4)` と略記すると、実 affine Pfaffian は順に

\[
\begin{array}{c|l}
A&n_0(n_9-w)\\
B&(n_9-w)\{n_0+n_5(-y+z^2/2)\}\\
C&-n_1(w-3n_9)^3/3\\
D&-n_6(w-2n_9)^2/2\\
E&-n_8\{x-yz+z^3/3+n_2+(3n_3/5)(w-n_9)\}\\
F&n_4(y-z^2/2+n_7).
\end{array}                                                   \tag{2}
\]

これにより表の rank四条件は必要十分である。
下記の未被覆四源 A,B,E,F は、四 raw Jacobi と実 Pfaffian を有限恒等式として検算した。
C,D は(1)への代入と次節の全支持帰着で扱う。

## 2. CとDは全て既存の冪零Jordan4成分に入る

五変数への斉次化後、全 modular 場は

\[
 N=x_2\partial_{x_1}+x_3\partial_{x_2}+h\partial_{x_3}.
\]

これは Jordan 型 (4,1) の nilpotent 線形場である。
C では ℓ=w−3n9h、D では ℓ=w−2n9h を、もう一つの超平面の式とする。
両者で n4=n5=n8=0 だから、全係数から

\[
 P^\sharp d\ell=-r\ell N,
 \qquad \Pi^\sharp d\ell=-(r+1/5)\ell N,
 \qquad \Pi=P-\tfrac15E\wedge N
\]

が成立する。従って ℓ を無限遠へ替えると、その affine modular 場は (5r+1)N となる。
ℓ の w 係数は1なので、ker N 上で ℓ=1 となる点を基点に選べる。
これにより定数部は零、線形部は rank三 nilpotent Jordan4 となる。
`MULTIPLE_HYPERPLANE_NILPOTENT_SUPPORT.md` の全支持定理が、全パラメータで適用できる。

さらに C の実因子は `2{h=0}+3{ℓ=0}`。
唯一の三重超平面 ℓ は `im N*` に入らないので、Jordan4 の s=1/5 成分に入る。
D の実因子は `3{h=0}+2{ℓ=0}` で、三重超平面 h は `im(N*)³` に入るから、
Jordan4 の s=1/15 成分に入る。
これは非線形 affine 自己同型を射影座標変換の代わりに用いる議論ではなく、
実際の一次超平面の変更と Euler 鉛直補正による全係数の帰着である。

## 3. EとFの全源は二つの22次元全成分になる

各源の全五自由パラメータを1にした実点を取り、全150二次 lift 係数と
全350射影 Jacobi を保った有限計算を行った。

| 源 | rank D の下界 | vertical | 鉛直・軌道・scaleの下界 | 全実源を加える下界 |
|---|---:|---:|---:|---:|
| E | 103 | 24 | 47 | 47 |
| F | 103 | 24 | 47 | 47 |

記録された非零 minor により下界は標数零でも成り立ち、D が全実方向を零化することは
有理係数で厳密に確認した。従って D の kernel 次元は高々47かつ少なくとも47である。
射影全接次元と実既約族の下界はどちらも22となり、点は滑らかで全成分は唯一。
各全源の PGL 飽和閉包はその22次元成分そのものとなる。
両代表の実軌道も22次元なので、それぞれの成分の稠密開軌道を与える。
数値の別実装による独立 QQ 再構成は別担当へ依頼中であり、producer の挟み撃ちと区別する。

### 実因子と相違

E の一般因子は `2H+C3`、F は `3H+Q3`。
E の三次式は x について一次で係数 h²、残部は h=0 で z³/3 となるので primitive。
Gauss の補題により C3 は既約である。
F の二次式 `hy−z²/2+n7h²` は rank三で既約である。
非被約五次因子という閉条件は両全源の閉包に成立し、SNC点を持つ MPS 成分とは異なる。
22以外の次元の成分とも異なる。
F は既存三重超平面成分の rank四二次式と異なり、他の22次元成分とは因子次数・重複度で区別される。

E と同じ `2H+C3` 型を持つ候補は、次の intrinsic な cofactor で区別する。
唯一の二重超平面を ℓ とし `Xℓ=Π♯dℓ/ℓ` とする。
この線形場の Jordan 型と、ℓ の `im Xℓ*` への所属は射影同型不変である。

| 成分の源 | cofactorのJordan型 | 二重超平面の一次式 |
|---|---|---|
| nilJordan4 c=0、s=1/3 | (4,1) | im Xℓ* に入らない |
| 本ノートE、nil(3,1)長い鎖top | (4,1) | im(Xℓ*)³ に入る |
| nil(3,1)短い鎖top、r=−1/3 | (3,2) | Jordan型だけで上二者と異なる |

上二者の代表はいずれも22次元の開軌道なので、射影非同型から成分の相違が従う。
短い鎖topの成分についても同じ開軌道の証明がある。
その源と proof は別担当の `MODULAR_NIL31_SHORT_SUPPORT.md` を用いる。

## 4. AとBの全支持は残す

同じ有限jobで A と B にもそれぞれ一つの全源代表を指定した。
結果は A が D rank下界98、軌道込み45、実源込み47、
B が D rank下界100、軌道込み46、実源込み47である。
よって実族の次元下界は22、周囲局所次元の上界はそれぞれ27と25。
この差を零とみなしたり、形式の有限次近似から全支持を被覆済みとしたりしない。
一般因子は A が4H+H'、B が2H+H'+Q3である。
新しい源は全式と全パラメータを保持しており、次は既知成分への実変形または真の支持判定が必要である。

## 有限検算の入口

`code/data/classification_nil31_long_sources.json` は未被覆だった A,B,E,F の全源と各一代表を保持する。
共通 runner `code/classification_modular_source_components.py` は全 Jacobi、全 source、Pfaffian、
全周囲微分を扱い、`CLASSIFICATION_NIL31_LONG_SOURCES_20261107.json` に保存した。
E_multiple_hyperplane の90秒・1536MiB契約 `jobs/nil31_long_sources.json` で実行、四源合計約1.66秒。
すでに幾何的に被覆した C,D に対する接空間の再計算は行っていない。

```sh
python3 scripts/classification.py gate research/classification/jobs/nil31_long_sources.json
ulimit -v 1572864
timeout 90s .venv/bin/python code/classification_modular_source_components.py \
  code/data/classification_nil31_long_sources.json /tmp/nil31-long-sources.json
```
