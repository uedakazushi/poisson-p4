# 残る三つの平面1+1+3全源：既知whole親はT36/T37/T38だけが残る

2026-11-08、PRIMARY E、`/root/log_manuscript`。
T21 marked-contact及びGa9 cuspidal-discriminantの二producer手証明を
別作者として独立読解して受理した。
併せて全既知126成分に対する下記の**whole-source親の必要還元**を手で得た。
この新全126比較の統合者による別読解は未了である。
新CAS、point、anchor、minor、normal jet又はGBは実行していない。
陽性の親、成分性又は全E被覆は依然未証明である。

## 1. 対象は三つの全六磁気sourceの閉包

base(h,x,y)、fiber(z,w)で二simple Hをh,xとし、三次を次の通り取る。

|source|C|
|---|---|
|node_flex_tangent_through_flex|h y²-x(x-h)²|
|cusp_flex_secant_through_flex|h y²-(x+h)³|
|cusp_cusp_tangent_cusp_secant|h²y-(h-x)³|

`TWO_LINES_PLANE_CUBIC_SUPPORT.md` §1の全source還元はfiber依存を残して
Saito基底へ戻したものであり、固定Dでは全GL2 anchor frameと全R2六係数を尽くす。
各canonical六磁気源のPGL5飽和閉包をS_iとする。S_iは既約である。

actualDの本質変数数は正確に三だから、そのconstant方向kernel Kの選択は
Gr(2,5)の6次元である。
quotient P2内で上表のmarked cubicと二線のPGL3 orbitは8次元。
実際connected stabilizerは二Hを個別に保つので、infinitesimalには

```
delta h=alpha h, delta x=beta x,
delta y=gamma y+d h+e x
```

と置ける。上の三Cそれぞれでdelta C=lambda Cの全係数を比べると、
alpha=beta=gamma、d=e=0、lambda=3alphaだけを得る。
これはscalar EulerなのでPGL3 stabilizerは零次元、従って有限である。
順序交換は有限二択で次元を変えない。

fixed projective Dの全Poisson源はGL2の4次元、magnetic六次元、
projective tensor scalarの一方向を引いて高々9次元。
この数え方ではfiber coordinate changesをもう一度足さない。
全h-Casimir代表と全GL2 Saito frameがfixed Dの全tensorを既に表示している。
よってdim S_i<=6+8+9=23。
既受理の全source微分rank48からdim S_i>=48-24-1=23なので、

```
dim S_i=23.
```

これらの微分階数を今回再計算したとはしない。
三つの全源の全D101/101/96、source48は既保存の全係数証明書を使う。

## 2. 既知23次元成分にwhole S_iは入らない

もしS_iが既知23次元成分Cに含まれれば、閉既約集合の次元からS_i=Cである。
以下はこの等号を除く証明であり、個別の特殊磁気点とCの交わりを除くものではない。

* MPS23二十九成分のそれぞれには、受理済みactual sourceの全edgeを零にした
  κ x0x1x2x3x4、κ非零のSNC essential5点がある。
  `ESSENTIAL_THREE_QUINTIC_FIRST_SEPARATE_REVIEW.md` §5で別Fraction確認済み。
  一方S_i全体はgradient rank<=3のclosed determinantal条件を満たすので等しくなれない。
* reduced二H+既約cubicの九既知23成分は、上表とは異なるmarked配置の稠密源を持つ。
  S_iと等しければ両稠密constructible像は同じ非空開で交わるはずだが、
  内在的な唯一の既約cubic、node/cusp、flex、二線の接触配置が異なるので不可能。
  全十三配置及び九既知源の相違は既受理のsupport noteを使う。
* 他の四十三の非MPS23成分は、generic actualDが上記のreduced1+1+3型と異なる。
  projective multiplicationのproper像、nonreduced条件、又はgeneric irreducibilityにより、
  reduced1+1+3型の稠密源S_iと等しくなれない。

次元23未満の十八成分もS_iを含めない。rank二の六成分はclosed rank条件により除外される。
これらは全sourceを一つの既知成分に入れる可能性の排除であり、境界の個別membershipの否定ではない。

## 3. T21全境界の二条件を全三源へ適用する独立確認

producer `TWO_LINES_PLANE_CUBIC_T21_MARKED_CONTACT_FRONTIER.md`
（SHA256 `58161610e7f82c2e305140539fdebd1602d16ae994bac4c51a5b61b387a29462`）を全読解した。
全period/edgeを残すT21 source及びproper ordered-two-factor flagから得る必要条件は

```
Z_h+X はx=0上でprojectively零、
Z_h+4X はh=0上でprojectively零、 X=P sharp dlog(x/h).
```

この元T21全五parameter Euler gaugeの係数は、既独立
`TWO_LINES_CUSP_TANGENT_SECANT_GLOBAL_PARENT_SEPARATE_REVIEW.md` で確認済み。
factor flagはproperで、target二Hがsimple/distinct、cubicがirredなので
極限のHはその正確な二因子から選ばれる。二条件はH交換でも交換される。

全targetではLiouville恒等式から
`Z_h=-P sharp dlog(xC)` である。
従って二場は `-P sharp dlog C` と `P sharp dlog(x³/C)` で、fiber方向だけを持つ。
hyperplane一般点でのprojective零はliteral零へ戻る。
x=0ではU_x,V_xがx倍、det(U,V)=unit*xCの単純因子条件からU_y,V_yの少なくとも一方が非零。
第一条件はC_y|x0=0を強制し、C|x0=c_h h³、c_h非零となる。
h=0ではraw4x4 matrixのPf=xCが一般点で非零なので、第二条件は
partial_x(x³/C)=partial_y(x³/C)=0を強制し、C|h0=c_x x³、c_x非零を与える。
ここでprojective rankとraw4x4のrankを同一視していない。

よってC=c_hh³+c_xx³+hx(d_hh+d_xx+d_yy)、d_y非零であり、
二Hは同じnodeの二接線となる。上表の三CはいずれもC|x0がh³の倍でない。
このため**全六磁気・零磁気を含む三源の全点がT21外**である。

## 4. Ga9 cuspidal24全境界の独立確認とscope明確化

producer `GA9_CUSPIDAL_MARKED_DISCRIMINANT_FRONTIER.md`
（SHA256 `fb2262422a15dfb2915a59a5d60a1d6cd5bd615cd3192ac9fe52bdecc898ad05`）を全読解した。
全sourceのactualDは hF_a、

```
F_a=(hW+6x²)²-144x³(h+ax).
```

exact essential3 locusではkernel Kとquotient P2が正則に復元される。
marked (H,F,p)のflagはproperで、target Hは二simple因子の片方に限る。
F|Hが非零fourth power、pで0/1jetが零、二次jetがH²のscalar倍という条件は全て閉。
最後の「二次cone」はこの0/1jet零も含む意味で用いる。
これを明示すれば、a0零を反転せず全flag座標で

```
F=a0 h²W²+hW B2(h,x)+C4(h,x)
```

となる。W^3 h等のlinear jetを残したままこの形を仮定するものではない。
binary G4=B2²-4a0C4はW→dW+Lでd²倍され、base linear変換にも三重根条件が不変。
G4=0を含む「l³m」のproper image coneは閉である。
親ではG4=576x³(h+ax)なので、全component boundaryにも必要条件となる。

上表の最初の二sourceで、親Hをflex線Zに選ぶとF|H=-X⁴だが、
pの二次coneはZ(X+Z)又はZ(X-Z)でH²の倍でない。
他方のHでは制限がZ²Y²又はZ²(Y-Z)(Y+Z)となりfourth powerでない。
第三sourceは唯一可能H=Yでa0=0、B2=hx、G4=h²x²となり、三重根が無い。
他方Hでは制限h³(W-h)がfourth powerでない。
従って**三つの全磁気sourceはGa9 cuspidal24の全体の外**である。
この結論はfixedsource slice、a0非零、quarticの一般点だけに狭めていない。

## 5. 全126の残りの高次元成分と正確な結論

参照prior126は公開commit `dab368089679508ae5c9016d66ba6accfe6d8fa9`、
catalogue blob `cd2e6c7bbd66c178777c85209253c5ff4fc10526`、
SHA256 `9e947a54b7cd1b7ef7cd48121080b54aaaaff82e65475828e3ec34acb873d5e6`。
全126行の必要fieldと理由をimmutable
`inputs/two_lines_plane_cubic_prior126_parent_projection_20261108.json` に保存した。

残る24次元のT02/T03/T16/T17、CPLUS_ext、rank3-doubleQ24、従属三H+Q24の七成分は、
generic Dの[1,1,1,2]又は[1,2,2]等のproper factorization条件に全体が入る。
targetには既約degree3因子があるため、そのpartitionには属し得ない。
25以上ではT01の五H条件、又は追加六成分の全nonreduced条件により同じく除外される。

T39/T40について、最初の二sourceは既受理D101の同じ点のtangent上界24を持つため、
whole S_iを25次元成分へ入れることは不可能。
第三sourceではI非零の全六磁気源を正則shear・実inverseで全lambdaへ戻した
既独立global134/45 graph判定が両方ideal(1)を与える。
従ってwhole S_iはT39/T40に入らない。
**I零の特殊磁気boundaryまでgraph陰性を延長しない**。

以上により、もし三S_iの一つが既知126成分の有限閉和に全体として入るなら、
その既約性から一つの成分に入り、そのwhole親は **T36、T37、T38のいずれか** でなければならない。
これは三つのいずれかに実際入ること、三者のいずれも親でないこと、又は新成分があることの証明ではない。
個別の特殊磁気点が別の23成分やT39/T40へ入る可能性も残す。

## 6. 停止済みuniform入口をpositiveへ読み替えない

保存されたT36 main uniform Lagrangian-plane restrictionは、非零scalarを除き

```
Q(B,C)*(3A Q(B,C)+25B³+5tC³),
Q=5tB²-tuBC+25C².
```

Qの二根が相異なれば、三次のp=[1:0:0]の二次coneもQなので、
二Hは同nodeの二接線である。Qがdoubleなら二Hも重複し、上表にはならない。
T38 mainはLQ²、T38 cell013 prime1の表示cubicはbinaryなのでC上三線へ分かれる。
これらは既保存の特定uniform入口のscope診断であって、whole T36/T38の陰性ではない。
全T37 B-curveのLQ²も既独立受理の限定的uniform入口である。

従って旧fixed-parent uniform normalcone、又は次数を増やすnormal jetの反復には戻らない。
次の必要な数学は、moving genuine periods/edgesを保つnonuniform actual confluence、
又は全known-parentへの別の構造橋である。三whole sourceの陽性帰属は未解決のまま保持する。
