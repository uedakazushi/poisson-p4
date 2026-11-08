# moving-simple m三の残二line＋五fixed：全vertex/Ore入口

2026-11-08、producer幾何。全三十七source certificateの残七labelだけを扱う。
全G閉鎖はまだ未証明。fullP2九を戻し、weight零Euler gaugeとpositive-weight gaugeを区別する。

## 全vertex選択

元coordinatesは(a,b,c,z,w)。全raw P0は三十七source reportのexact式を使う。
二whole line DEG1_HA_CA_INTERSECTION/HB_CA_INTERSECTION、
fixed DEG1_TRIPLE_EXACT_B及びDEG2_CB2_1はvertex cで元P0がc非依存。

DEG1_HA_CA_OTHERのt=−8/3では唯一のc依存項は−3cz∂c∧∂w。
weight零Euler Lw=−3zを引いてc依存を消す。
y=(a,b,z,w)の全Q0は

`Q0=(6az∂a+3bz∂b+3z²∂z)∧∂w`。

DEG1_HB_CA_OTHERのt=4/3ではLw=3zで同様にc依存を消し、

`Q0=(−3az∂a−6bz∂b−3z²∂z)∧∂w`。

いずれもw非依存の定数∂w因子を持つ。
ただしhigher Qまでこの性質が延長するとはまだ主張しない。

DEG2_CA2_0ではvertex bを使う。
L_c=2aをEuler gaugeとして引けば全b依存が消え、P0のb-row全零。
y=(a,c,z,w)上で

`Q0ac=−3a², Q0aw=6cz, Q0cw=−6aw`、他のbracketは零。

その二CasimirはzとI=a²w+c²zである。
weight零W0は零。元leadingはz²Iである。
raw liftのbase Euler gaugeを使ったので、rank零planeのFrobenius補題を流用していない。

## 全positive-weight source

各sourceのfullP1必要kernel及びFULLP2九に対し、選んだvertex vを使い
L_i=−∂vΠ^{i,v}、L_v=0を引く。
その後残normal-quadratic Euler冗長性二でWzのz²、Wwのw²を消す。
全sourceがv非依存に入り、Q(y)+∂v∧W(y)となるかを全係数で検証する。
入らないsourceがあればそのexact obstructionを残し、wholeOreを仮定しない。

全Q係数空間のrank/basisと全W十八cover、元sourceからQ/Wへの全linear mapを保存する。
Q0非零のsourceではQ/W十八を独立に取ることは一般にsupersetであり、
actualD3の元line条件や他のsource相関を使う際には元linear mapを必ず戻す。

元projective Jacobiのvertex次数からstrict[Q,Q]零とE_y∧[Q,W]零が必要。
次の一有限batchは全Qspace/strictJac及びactualD3 coefficient equationsを作る。
旧negative constant-rowspan factor criterionをもう一度回さず、
全strict-support/actualD/normalizerの構造へ進む。

generic triple、m二、reduced R、他の未独立pointwise Eは別義務。
