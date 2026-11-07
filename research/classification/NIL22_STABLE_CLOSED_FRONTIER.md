# stable 二源に対する二つの既知24次元閉包の除外

2026-11-07。全六磁気係数を保持する nil22-image の stable_line_conic と
stable_three_lines を対象とする。証明済みの全十一源表示と
`NIL22_STABLE_PROJECTIVE_ACTION.md` の因子配置を入力とする。
本稿は CPLUS_EXT と dependent-3H+Q の既知24次元成分だけを除外する。
別成分への帰属、新成分性、E 全体の完了を結論しない。

## 1. CPLUS_EXT の閉じた共通接線条件

既知 CPLUS_EXT の一般因子は `H Q1 Q2` であり、平面商で H は両 conic の
同一点での共通接線である。従って H 上の二つの二次式の制限は、
同じ一次式の平方の倍数となる。

H,Q1,Q2 の射影パラメータ空間で、両制限の rank≤1 と、二制限の係数行列の
rank≤1 を課す。零制限も許すので、これは閉条件である。
因子の乗法写像は proper なので、この条件を満たす因子分解を持つ五次式全体も
閉集合となる。従ってこの条件は CPLUS_EXT の全 rank4 閉包で成立する。

stable_line_conic の五次式は `h² L Q`、L は Q の接線、h は割線で、
接点は h の外にある。外側 H の選択は h 又は L に限られる。

- H=h なら二次因子は Q と hL。Q の制限は二つの異なる一次式の積なので rank2。
- H=L なら二次因子は Q と h²。両制限は平方だが零点が異なり、比例しない。

stable_three_lines は `h² L1 L2 L3` で、四直線のどの相異なる三本も独立である。
H=h なら残る一つの二次因子が相異なる二 Li の積になり、H 上で rank2。
H=Li なら、残る二次因子対は `(h²,Lj Lk)` 又は `(h Lj,h Lk)` である。
どちらも少なくとも一方の H 上の制限が rank2 になる。

従って両源の全点は、この CPLUS_EXT の必要な閉因子条件を満たさない。
磁気係数は因子式を変えないので、その例外もない。

## 2. 従属三超平面の二重化に関する regular collision 補題

実 rank4 Poisson 族の一般退化因子に、相異なる三超平面があり、
その一次式が二次元を張るとする。特別点でそれらが `h,h,L` へなり、
h,L は独立と仮定する。族の base を有限拡大して因子をラベルし、
独立な対 h_t,L_t を regular な線形座標変更で h,L へ固定できる。
残る第三因子は

`h'_t=α(t)h+β(t)L`, `α(0)≠0`, `β(0)=0`

となる。一般点では三因子が異なるので β は恒等零ではない。
h-Casimir 代表を取る操作は固定座標での線形 Euler 商操作なので regular な族を保つ。
affine h=1、x=L と置く。

L と h'_t がともに Poisson 因子であることから、一般 t で全座標 z_j について

`{x,z_j}_t=x(α+βx)c_j(t)`

となる。実際、両因子は互いに素で、左辺の次数は高々二なので、商は定数である。
x の係数は `α c_j`。左辺の regular 性と α の unit 性から c_j も regular となる。
従って極限では

`{x,z_j}_0=α(0)c_j(0)x`。

すなわち h=1 における対数 Hamiltonian `P#dL/L` は translation であり、
その線形部分は零になる。これは正規化した座標の選択によらない。
cofactor を先に反転して非正則な極限を採ったのではない。

この補題は従属三 H と Q の全源閉包の rank4 点に適用できる。
Q が極限で可約となることを許し、上の三 H の選択可能性を全て保持する。

## 3. 全因子選択と全六磁気源への適用

stable_line_conic で三 H を選ぶなら `h,h,L` しかない。
stable_three_lines で相異なる三 H を選べば一次式の rank は3であり、
従属という必要閉条件に反する。残りの全選択は `h,h,Li` である。

一方、両 stable 源の全ての単純一次因子 L について、
`P#dL/L` は磁気係数に依らない affine 場であり、その線形部分の rank は厳密に2。
全係数、平方零性、rank を下記 checker に保存した。
従って §2 の translation 必要条件に反する。
以上により両完全六磁気源は、既知 dependent-3H+Q の24次元成分の外にある。

## 4. 有限再計算と射程

`code/classification_nil22_stable_closed_frontier.py` は両源の全 factor assignment を列挙し、
CPLUS の制限 rank/比例 rank と、全単純因子の対数 Hamiltonian の全線形係数を再構成する。
`jobs/nil22_stable_closed_frontier.json` の有限 gate の後に実行する。
出力は `reports/CLASSIFICATION_NIL22_STABLE_CLOSED_FRONTIER_20261107.json`。
§2 の regular collision 補題は上記の手証明を用い、有限標本や数値極限で代用しない。

元の stable 源とその PGL 飽和は既約であるが、本稿はその成分性を証明しない。
二つの既知閉包以外への非所属も主張しない。
