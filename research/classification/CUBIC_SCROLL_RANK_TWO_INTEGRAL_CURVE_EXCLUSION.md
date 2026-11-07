# smooth cubic scroll・rank 二 restriction・既約被約零曲線の排除

2026-11-07。主義務 G。root の幾何提案を別作者が再導出し、正規束の二つの制限と
全 normal-cone 分岐を直接確認した。新しい全計算の成功から推測する命題ではない。

**命題。** `S=S(1,2)≅F1⊂P4`、`H=C0+2f` とする。rank 四 Poisson tensor の実退化因子
`D` が既約被約 quintic、`S⊂Sing D`、`pi|S` の一般 rank が二であるとする。
その零因子 `R=div(pi|S)∈|-K_S|` が integral reduced なら、この source は空である。
`R` の smoothness、横断 ordinary node、vertical 成分の不在を仮定しない。
可約・非被約 `R` と cubic cone はこの命題の補集合であり、全 G の解決ではない。

## 1. intrinsic abelian conormal と parallel quadratic

nonplane `S` では secant argument により `D` の一般 normal multiplicity は二である。
`S` の一般 smooth point において、接方向に非退化な rank 二部分を取り、normal
座標を `u,v` とする。normal bracket の一次項が非零なら、全 Pfaffian の一次項は
接方向の非零 Poisson coefficient とその一次項の積となる。`D∈I_S^(2)` に反する。
従って intrinsic conormal bracket は abelian、すなわち `{I_S,I_S}⊂I_S²` である。
これは rank 零の場合の Frobenius-Lie 議論をそのまま転用したものではない。

`U=S\R` では conormal abelian Lie algebroid の Hamiltonian action が `N*` に
flat connection を与える。Jacobi によりこの作用の curvature は零であり、dual の
正規束 `N=N_(S/P4)` も flat である。`sigma=pi|S` とし、実五次式の normal quadratic を

`G∈H0(S,Sym²N*⊗O_S(5H))`

と書く。adjunction は `det N=K_S+5H` を与えるので、`G/sigma` は
`Sym²N*⊗det N` の parallel section である。この主張は Hamiltonian flows が全 `pi` と
`sigma` を保つことを normal quadratic へ制限しても得られる。
従って `G` の rank は連結な `U` で一定であり、その零・rank 降下 locus は `R` 上だけである。
また

`disc(G)=c sigma²`, `c∈C`

が全 `S` 上の section 恒等式として成立する。実際、左辺を `sigma²` で割ったものは
parallel な有理関数であり、rank 二 anchor が `T_U` 全体を張るため定数である。
一般 normal multiplicity 二より `G` 自体は零でない。

## 2. rank 一の normal quadratic

`c=0` とする。generic rank 一の quadratic を height 一で primitive に分解し、
line quotient `ell:N→A` に対して `G=gamma ell²` と書く。
`S` は smooth なので飽和した rank 一因子は line bundle に延長し、primitive 性により
`gamma` の divisorial pole はない。`U` 上で nonzero parallel なので
`V=div(gamma)` は `R` に支持される。`R` が integral reduced だから `V=kR`、`k≥0`。
従って

`2A+kR=5H`。

`R=-K_S=2C0+3f`、`H=C0+2f` を代入すると、右辺から `kR` を引いた `C0` 係数は
`5-2k` である。これは奇数なので、整数 line bundle `A` の二倍にならない。矛盾する。

## 3. rank 二と全 vertical 分岐

`c!=0` なら `disc(G)` は `C(S)` 上で平方なので、二つの generic root は rational に分離する。
それぞれを primitive な line quotient `ell_i:N→A_i` に延ばし、
`G=gamma ell_1 ell_2` とする。`delta=ell_1 wedge ell_2` は
`A1+A2-detN` の正則 section である。判別式恒等式から

`V+div(delta)=R`, `V=div(gamma)`, `V≥0`, `div(delta)≥0`。

ここで `R` が被約既約なので、全場合は `V=0` 又は `V=R` の二つだけである。
vertical case を除外して始めてはいない。

### 3.1 `V=0`

`Y=Bl_S(P4)`、`E=P(N)` とする。tangency と intrinsic conormal abelian により `pi` は
この blowup に正則に持ち上がる。実反標準因子は `D'+E` であり、inverse closed 2-form の
`E` residue は closed rational 1-form `alpha` である。その pole は二つの horizontal prime
`Gamma1,Gamma2⊂E` だけで、generic fiber 上の residue は非零定数 `lambda,-lambda`。
これは `c!=0` と generic normal quadratic の相異なる二根から直接分かる。
`V=0` なので余分な vertical pole はない。

全点の SNC を仮定せず、必要な数値 Gysin 関係だけを得る。`xi=O_E(1)` とし、
十分大きな `n` に対する ample divisor `A_n=xi+n p*H` を使う。
`A_n²` と `A_n(A_n+p*L)` の一般 complete-intersection 曲線を選ぶ。
`L=H,H+f` はともに very ample として取れる。これらの一般曲線は二 prime の交わりや
各 prime の singular locus 等の codimension 二 bad set を避ける。曲線上の residue 定理より
`lambda(Gamma1-Gamma2)·C=0`。
`[Gamma1]-[Gamma2]=p*(A1-A2)`、`p_*xi=1` として二種の曲線を引くと

`(A1-A2)·L=0`, `L=H,H+f`。

よって `A1≡A2`。`Pic(F1)=Z C0⊕Z f` なので line bundle 自体が同じである。
一方 `V=0` から `A1+A2=5H` であり、`C0` 係数五は二倍に分けられない。矛盾する。

### 3.2 `V=R`

この場合 `delta` は零点を持たず、`N≅A1⊕A2` が全 `S` 上で成立する。
次節の elementary normal-bundle 計算によりこれは不可能である。

## 4. 正規束は分解しない：二本の line と c2 の直接計算

adjunction と自交差から

`c1(N)=3C0+7f`, `c2(N)=deg(S)²=9`。

任意の ruling line `f` について normal sequence は

`0→O_(P1)→O_(P1)(1)^3→N|f→0`。

左の nowhere-zero map を与える三つの linear form は二次元を張る。よって Euler sequence と
余る一 summand から `N|f≅O(1)⊕O(2)`。

negative directrix `C0` は parametrization
`(su,tu,s²v,stv,t²v)` の `v=0` という line である。
normal sequence の左写像は

`0→O(-1) --(s²,st,t²)--> O(1)^3→N|C0→0`。

三つの quadratic の二つの linear syzygy `(t,-s,0),(0,t,-s)` により、その quotient は
`O(2)⊕O(2)` である。特に単に degree 四から splitting を推測したのではない。

もし `N=A⊕B` なら、ruling への制限から順序を取り直して `A·f=1`、`B·f=2`。
directrix への制限から `A·C0=B·C0=2`。従って唯一の候補は

`A=C0+3f`, `B=2C0+4f`。

しかし `A·B=2C0²+10C0·f=8` であり、`c2(N)=9` と矛盾する。
よって `N` は分解せず、§3.2 も排除された。これで全 normal rank と全 vertical 分岐が閉じ、命題が従う。

## 5. 完全線形源と残る frontier

別の producer 有限計算は、全150係数の tangency と、正確な intrinsic bracket

`H_ij=x0{qi,qj}-2qj{qi,x0}+2qi{qj,x0}`

の全 first partial を dominant parametrization に pullback した。
`H_ij` を raw `{qi,qj}` に置換してはいない。全 tangency rank は87、結合線形系は
`1066 x 150`、rank105、kernel45。全 Euler24 を確認した canonical 源は完全な21次元である。
全350 Jacobi は独立39二次式、実Dの係数のその二次 rowspace への包含は false。
一回の弱い leading-monomial 上界は固定源15、全射影像32であり、G は閉じない。

`Aut(S)` は connected、`S\C0` に推移的である。`p=(0,0,1,0,0)`、chart `x2=1` では
接座標 `k=x0,r=x3` における restriction の係数は canonical source の `t2` 一個となる。
全 rank 二 component は actual Aut により `t2!=0` patch に入り、scalar で `t2=1` と出来る。
元350式の rows `80;26,97;81;57;27,99,117;59` は九変数を一定係数で解く。
対応する非零 determinant は `-1,-15/2,-3,8,-240,2`。
全 source を戻した残十一変数への完全350式の弱い上界は normalized dimension10、全像28であり、
これも非優越証明ではない。目標は normalized dimension **六** 以下である。

二回の弱い dimension 判定は未閉として保存し、GB/jet/row-order 反復をしない。
数学的な次の方向は、上で証明した integral reduced `R` の排除を使い、可約・非被約の
anti-canonical `R` と実Dの既約性を全源で扱うことである。
補助的に rows73,80 は正確に

`t18=(2t0-t17)(2t0-2t17+3t19)/3`,

`(t0-t17)(2t0-t17)(2t0-2t17+3t19)=0`

を与え、全源の三つの linear branch が得られる。しかしその三 branch と `R` の因子型の
完全対応、各 full rank 四源の非優越はまだ証明していない。

有限計算は producer 証拠であり、この幾何命題の独立証明と区別する。

- runner: `code/classification_cubic_surface_tangent.py`。旧七 FunctionDef AST を保つ isolated open mode。
- 全線形源: `reports/CLASSIFICATION_CUBIC_SURFACE_TANGENT_SMOOTH_SCROLL_12_20261107.json`。
- 全 actual patch: `reports/CLASSIFICATION_CUBIC_SCROLL_RANK_TWO_OPEN_20261107.json`。
- 有限 job: `jobs/cubic_surface_tangent_preflight.json`, `jobs/cubic_scroll_rank_two_open.json`。
- 元の preflight job bytes は `inputs/cubic_surface_tangent_preflight_job_snapshot_20261107.json` に不変保存。

preflight 2.182秒、RSS56,904KiB。全rank二 patch 18.431秒、RSS78,612KiB。
両 gate は PASS だが、両 dimension goal は false として保存した。新しい成分数、全 G/E の完了、
cubic cone の排除又は全 `R` の空性は主張しない。
