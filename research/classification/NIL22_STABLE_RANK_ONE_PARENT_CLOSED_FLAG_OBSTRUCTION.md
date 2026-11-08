# Stable nil22 二全源：rank一affine modular親の閉flag障害

2026-11-08、PRIMARY E_multiple_hyperplane。純手補題。新operator／point／profile／jetは零。
本稿は新130の非定数g long source、および新三共鳴・binary-cubic長鎖候補を、
両stable全源の親候補から一括して除く。whole E、T36/T37非所属又は新成分性は結論しない。
独立読解前のproducer記録である。

## 1. marked重複Hに対する閉量

rank四のprojective Poisson tensor πと、実退化因子Dに対しH²|Dとなるmarked超平面Hを考える。
Dの因子はPoisson divisorなので、Hを無限遠にしたquadratic affine代表が存在する。
そのH-Casimir homogeneous lift P_Hは一意である。
差がE wedge YでH行零ならH Y-Y(H)E=0からそのwedgeは零となるためである。
Hのある係数が非零の各chartでHを第一座標に取ると、
Euler gaugeは線形係数抽出と正確な一次式除算で正則に作れる。
従ってcolumn divergenceは正則なaffine vector field

```
Z_H=c_H+A_H*x
```

を与える。A_Hは対応する斉次場のkerHへの制限である。
固定Hを保つprojective変換はaffine変換なので、線形部はGL4共役で移り、
非零tensor尺度はA_H全体をscaleする。よってrank A_Hは座標・尺度に依らない。
斉次場の全rank、Jordan型又はaffine定数部をこの線形部と混同しない。

各H-chartでrank A_H<=1は全2×2 minorsの消滅という閉条件である。
因子flag D=H²*C3をprojective係数積P(H)×P(C3)で取り、
Poisson接性と上のrank条件をその上で課す。
rank四openではD!=0なのでtensorから[D]への射は正則。
projective multiplicationのpropernessから、flagを忘れた像もこのrank四openで閉である。

従って、全点にrank A_H<=1のmarked repeated Hを持つsourceのPGL飽和閉包は、
任意のrank四極限点でも**何らかの**marked repeated Hで同rank条件を持つ。
marked Hがtargetのある任意のsimple factorへ移ったと仮定してよいわけではない。
H²|Dというmultiplicity条件もproper flagに残る。

## 2. 全rank一long sourceに適用する

`RANK_ONE_LONG_NONCONSTANT_COMPONENT_ENTRANCE.md` とその独立reviewの全源は

```
X=∂x+x∂y,
N=h∂x+x∂y
```

をcanonical modular場に持つ。固定hのaffine線形部はx∂yでrank一、
斉次NはJordan型(3,1,1)、斉次rank二である。
この全sourceはDにh²を持つ。全parameter、全magnetic lower係数、実affine変換、
全PGL飽和と非零tensor尺度を戻しても§1のrank<=1条件を満たすmarked Hを持つ。
従って、新130の全閉成分で同条件が保持される。
この結論は一つのgeneric pointからの外挿ではなく、その全sourceの恒等的modular fieldを使う。

同じXを持つ `RANK_ONE_LONG_RESONANT_AND_BINARY_CUBIC_ENTRANCE.md` の
μ=1、μ=1/3、g零scalar-U binary-cubicの三whole sourceにも、そのまま同じ補題を適用できる。
これら三像のambient成分性を仮定する必要はない。
全sourceの閉包に成立する必要条件だからである。
短鎖rank一translationなど、marked Hに対しrank A_H<=1とh²因子を保つ他源も同じ結論を持つ。

## 3. 両stableの唯一flagと全磁気modular field

`NIL22_STABLE_PROJECTIVE_ACTION.md` の両全六magnetic sourceは、順序(h,q1,q2,x,y)で

```
P_h*=0, {qi,pj}=Mij, {x,y}=0,
{q1,q2}=R+h(V2*q1-V1*q2), p=(x,y).
```

全column divergenceを取ると両源とも厳密に

```
Z_h=x∂q1+y∂q2.                                           (1)
```

Rはqに依らず、δZに寄与しない。
line-conicではdivM一行目=-x+3/2、二行目=-y-3/2、
V=(3/2,-3/2)だから定数項は相殺する。
three-linesではdivM=(-x,-y)、V=0で同じ(1)となる。
したがってaffine線形部は全磁気係数でrank二、平方零、定数部零である。

実Dは全磁気係数から独立で

```
line-conic: h²*(x-y-12h)*(36h²+3h(x-y)+2xy),
three-lines: h²*(x+3h)*(y-3h)*(x-y-3h).
```

前者のquadricはrank三でC上既約、後者の三simple lineとhは相異なる。
両方とも唯一のrepeated hyperplaneはhで、重複度二。
§1のproperflag極限のHはしたがって**hだけ**である。
他のsimple lineへの交換、quadric因子との交換、動くHの取り直しではrankを回避できない。
しかしhで(1)のrankは二なので、rank<=1と矛盾する。

## 4. 結論と次の正の義務

両stableの全六磁気source、その全PGL軌道は、
新130および上記三つのfixed-long-X候補の全PGL source閉包の外にある。
新三候補の成分計数が未確定でも、このwhole source障害は有効である。
磁気不変量を反転せず、零磁気、zero-residue境界も含む。

単にJordan closure (3,1,1)→(2,2,1)が不可能と言ってはいけない。
通常のGL5 nilpotent orbit closureではその型の変化自体は可能である。
今回不可能なのは、multiplicityにより固有のHを保ったparabolic/affine問題での
線形部rank一からrank二への増加である。

本補題は新130・新三候補という有限の親候補に必要な一つの閉flag検査であり、
追加のnegative profile系列を開始しない。
R4_nil31の既存高次元親との因子閉条件、既T36/T37へのactual confluence、
未知親又は全成分性は別の数学義務である。
