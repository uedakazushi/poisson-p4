# λ=0・joint半単純rank一：全H＋二quadric像は次元24以下

2026-11-07。課題G及びE_reduced_nonnormalの入口。全因子型は
`D=hq₁q₂`、qᵢは相異なる既約二次式、`Pf P=κq₁q₂`, κ≠0。
二normal場Xᵢは `TWO_QUADRICS_SUPPORT.md` 第2節で定義し、ここではλ=0、
joint半単純部のrankが一とする。結論は**未知W優越成分の排除**である。
各源が既知成分のどれに入るかというEの全認識は主張しない。

## 1. 追加nil方向も保持した五つの完全supersource

二場の共通半単純生成元をSとする。qᵢはSで不変、trace S=0、hは固定である。
rank三以上のweight零二次形式があるため、五weightは `0,a,−a,b,−b`。
これはrank二の場合と同じ反対weight対・traceの議論で従う。
全weightが相異なれば可換nilpartは零となり、X₁,X₂はSの倍数だけになって矛盾する。
従って、非零scaleと座標置換を除き、Sのaffine weightは

`(1,1,−1,−1)` 又は `(1,−1,0,0)`。                        (1)

`Z=X₁+X₂` のJordan分解を考える。そのnilpartが非零ならそれをNに取る。
零なら、Xᵢのうち非零のnilpartをNに取る（なければ両Xᵢは従属）。
Jordan分解は有限次元polynomial表現にも適合するので、いずれの場合も

`[S,N]=0`, `Lie_S P=Lie_N P=0`, `div P=sS+tN`             (2)

となる。Xᵢのnilpartが互いに独立であっても構わない。
**Xᵢ自身がspan(S,N)に入るとは仮定しない。**

第一weight型のcentralizer内の非零nilは、二つの二次元weight空間内のJordan2の
片方だけ又は両方。第二weight型では、hを含む三次元零weight空間上にnilが作用する。
hを固定したJordan2ではhがimage内か外か、Jordan3ではhがその唯一kernel直線となる。
従って全て次の五型へ帰着する。

| 型 | S | N | 全200×92線形系のrank | 自由係数 |
|---|---|---|---|---|
| A1 | x∂x+y∂y−z∂z−w∂w | y∂x | 81 | 11 |
| A2 | 同左 | y∂x+w∂z | 83 | 9 |
| B1 | x∂x−y∂y | w∂z | 80 | 12 |
| B2 | 同左 | ∂z | 79 | 13 |
| B3 | 同左 | w∂z+∂w | 83 | 9 |

Pの全90係数にs,tの二変数を加え、(2)の20+180行を全て解いた。
全92係数の完全kernel、非零rank minor、全四Jacobi、実Pfaffianを保存した。
これは必要条件のsupersourceであり、残るJacobi/因子条件を勝手に消していない。
以下のn番号は `CLASSIFICATION_TWO_QUADRICS_MIXED_NORMAL_SYMMETRIES_20261107.json`
各型の完全源の番号である。

## 2. A二型の全支持

A1の**全**Pfaffianはyで割れる。従って二既約quadricの積にはならず、全型が除外される。

A2では `s=n₇,t=n₈`。実Jacobi係数 `−(2n₅+n₇)²` から被約支持で
`n₅=−n₇/2`。s≠0なら、実係数 `−n₃(3n₅+n₇)` からn₃=0、
次いで実係数 `−n₀n₅−n₀n₇+n₂n₃−n₃n₄+n₃n₈` からn₀=0。
この三代入後の全Pfaffianはywで割れ、対象型ではない。

s=0なら `div P=n₈N`。n₈=0で非退化Pなら `P#d(Pf P)=(Pf P)div P`
からPf Pが定数となって対象型に矛盾する。従ってn₈≠0、全tensorのscaleで
canonical nil22-imageへ入る。
`TWO_QUADRICS_NIL22_BASE_SUPPORT.md` の全源定理により、被約H+二quadricは必ず平面七配置。
その全PGL像の次元は高々24である。ここでも特定の平面源の成分認識は使わない。

## 3. B1の必要支持は六parameter

S,Nで不変な各quadricは

`qᵢ=aᵢxy+bᵢ(h,w)`、bᵢはbinary quadratic。

rank三以上からaᵢ≠0なので、各因子をscaleしてaᵢ=1とする。
Sの非零weight空間x,yは一次元なので、追加nilpartの有無に関係なく
`Xᵢ(x)=αᵢx`, `Xᵢ(y)=−αᵢy`。
全源のx normal条件から、両αᵢはn₁₀−2n₅である。
その和がs=n₁₀だからn₁₀=4n₅、αᵢ=2n₅。
joint半単純rank一からn₅≠0である。
Pfのu²係数（u=xy）は−2n₅n₈。両qᵢのu係数が一なのでn₈≠0。

元のJacobi係数を順に用いると、

`n₂=0`, `n₇=−n₁n₈/n₅`, `n₆=−n₀n₈/(2n₅)`, `n₀n₁=0`。

さらに定数係数は `n₀(n₁₁+2n₉)−2n₁²n₈/n₅=0`。
n₁≠0ならn₀=0となって矛盾するため、n₁=n₇=0。
各bᵢのnormal条件は

`−n₀+n₅w bᵢ'(w)=2n₅ bᵢ(w)`、従って
`bᵢ=Aᵢw²−n₀/(2n₅)`。

全Pfのw⁴係数は零なのでA₁A₂=0。
n₀=0なら一方のquadricはxyとなって可約であるため、n₀≠0。
従って最後にn₁₁=−2n₉を得る。
全対象支持は、六自由係数

`n₀,n₃,n₄,n₅,n₈,n₉` （n₀n₅n₈≠0）                   (3)

の完全源内に入る。全残余Jacobiはこの源で恒等零。
因子条件をさらに課せば次元は減るだけである。

## 4. B2の必要支持は二つの七parameter源

再びqᵢ=xy+bᵢ(h,w)であり、両αᵢ=n₁₁−n₆。
発散からn₁₁=2n₆、αᵢ=n₆≠0。
Pfのu²係数は−n₆n₉なのでn₉≠0。
実Jacobi係数から

`n₂=n₁₀n₆/n₉`, `n₁n₉=2n₁₀n₅`, `4n₁₀n₂=0`。

被約支持でn₁₀=0、従ってn₁=n₂=0。
残る全Jacobiは正確に

`n₀(n₁₂+n₈)=0`, `n₀n₉−n₅n₈+n₆n₇=0`。

n₀≠0ではn₁₂=−n₈とn₇=(n₅n₈−n₀n₉)/n₆を使い、
n₀=0ではn₇=n₅n₈/n₆を使う。
前者の自由係数は `(n₀,n₃,n₄,n₅,n₆,n₈,n₉)`、
後者は `(n₃,n₄,n₅,n₆,n₈,n₉,n₁₂)`。
どちらも七parameterで、n₆n₉≠0以外のminor補集合を捨てない。

## 5. B3の必要支持は二つの五parameter源

`v=w²−2hz` とすると、全S,N不変二次式は `a xy+bv+ch²`。
もし両qᵢのxy係数が非零なら、normalize後の両αᵢはs=n₇となり、
和の条件からn₇=0、両αᵢ=0となってjoint半単純rank一に矛盾する。
両係数が零でもαᵢ=0となる。
従ってちょうど一方だけが非零で、因子を交換・scaleして

`q₁=xy+Av+Ch²`, `q₂=v+Dh²`。

q₂のx normal条件は `−2n₂/(v+D)` が定数であることからn₂=0、α₂=0。
α₁=n₇≠0。実Jacobi係数−n₅n₇からn₅=0。
この時点で全Pfaffianは正確に

`(n₀+n₁v−n₇xy)(n₄+n₆v)`。

q₂が非定数なのでn₆≠0。全残余Jacobiは

`n₁(4n₆+n₈)=0`, `n₀(2n₆+n₈)+2n₁n₄=0`。

n₁≠0ならn₈=−4n₆、n₀=n₁n₄/n₆で、
自由係数は `(n₁,n₃,n₄,n₆,n₇)`。
n₁=0ならq₁の既約性からn₀≠0、n₈=−2n₆で、
自由係数は `(n₀,n₃,n₄,n₆,n₇)`。
両枝とも五parameterであり、全Jacobiは恒等零である。

## 6. normalizerの自由作用から全PGL像を評価する

B各型では、次のaffine群Kがhを固定し、Sを固定、Nの直線を保つ。
x,yの独立な二つの非零scaleに加え、(z,w)へ次を作用させる。

| 型 | (z',w') | Kの次元 |
|---|---|---|
| B1 | (az+bw+e,dw)、ad≠0 | 2+4=6 |
| B2 | (az+bw+e,dw+f)、ad≠0 | 2+5=7 |
| B3 | (a²z+ab w+c,aw+b)、a≠0 | 2+3=5 |

従ってKは各完全supersourceと実対象subset Tを保つ。
対象Tがnormal場の選択に由来する場合にも、S固定・N直線保存のためこの性質は保たれる。
`PGL₅×P(T) → P(全Poisson空間)` のfiberには
`(g,P)↦(gk⁻¹,kP)` というKの自由作用がある。
第一成分への作用が自由なので、Pのstabilizer次元を別に足す必要はない。
Tのaffine次元上界dはB1/B2/B3で6/7/5、Kの次元も同じである。
従って各全PGL像の次元は

`24+(d−1)−dim K ≤23`。                                  (4)

A2の残余像だけは既済の平面全支持から高々24。
以上でjoint半単純rank一の全像は有限個の次元24以下の構成可能像に含まれる。
その閉包も次元24以下なので、未知W優越成分（必要次元25–27）の一般元にはなれない。
既知成分への全帰属又は全源の局所唯一性はこの次元論からは従わない。

## 7. 算術状態と未完了部分

producer `code/classification_two_quadrics.py` の
`--mixed-normal-symmetries` と `--mixed-normal-support`、契約
`jobs/two_quadrics_mixed_normal_{symmetries,support}.json`、reports
`CLASSIFICATION_TWO_QUADRICS_MIXED_NORMAL_{SYMMETRIES,SUPPORT}_20261107.json`
は、全200×92源及び上記全stratumの元Jacobi/Pfaffianを検証した。
独立 `code/classification_two_quadrics_review.py --mixed-normal-support` は
producer/CASをimportせず、全五型の全200×92行、完全kernel、非零rank minor、
全四Jacobi・実Pfaffian及び全枝のLaurent代入を再構成してPASS（0.252秒）。
証明書は `CLASSIFICATION_TWO_QUADRICS_MIXED_NORMAL_SUPPORT_REVIEW_20261107.json`、
実行前gateは `jobs/two_quadrics_mixed_normal_support_review.json`。
有限五型の幾何完全性、normal条件からの必要strata、自由K作用による次元論は
本稿に示した別の証明入力であり、算術reviewだけが証明したとは扱わない。

λ=0のjoint半単純rank零（二normal場の全nilpotent配置）は未処理である。
rank一でのE全成分認識、及び平面tan-simple/chordのE帰属も残る。
