# 全ての Pκ 上で交叉は滑らかな R19 の通常二重構造である

2026年11月3日。これは元のR八係数でt*q5*q7≠0の部分と、その全PGL5飽和に対する定理である。
AF=0やt=0の全境界点へ広げない。

## 定理

全てのκ∈Cについて、標準形Pκの近傍ではCA,CB,Rはいずれも滑らかで、実交叉D=CA∩CBの完成局所環は

    C[[u1,…,u19,epsilon]]/(epsilon²)

に同型である。
従ってこの部分では交叉はCohen–Macaulayな通常二重構造で、associated primeはRだけであり、embedded primesはない。
これは全Poisson環の非CM性とは異なる。Pκを通る全Poisson環はLPB39とCA/CB28の存在によって非CMである。

## 1. 全κで消えない四つの実minor

既存の全canonical150係数微分行列について、κ=0で選んだminorをκを記号のまま計算すると、次の非零定数となる。

|行列|サイズ|全κでの正確な行列式|
|---|---:|---|
|CAの全source|29|3623878656/390625|
|CBの全source|29|8153726976/78125|
|両tangentの和|37|179381993472/390625|
|Rのpaired-flag source|20|9663676416/390625|

一般κでの上界は既存のQQ(κ)上の全行列計算とsource次元から得るので、
全κでrank29/29/37/20が厳密に一定である。
有限個の標本を補間したものではない。

## 2. 全κの全旗fiberとproperなjoint source

Q=h²+κhz+z²は常に非零である。
既存のcofactorによる全旗分類はκに制限を置かず、CA旗とCB旗をそれぞれ一意に戻す。
Pκの非N性も、唯一のself-cofactor line xとNの第二flag条件の不一致から全κで従う。
両proper模型のsource微分は全κでprojective rank28なので、
全fiberは長さ1となり、有限化とNakayamaで両成分の全近傍と同型となる。
従ってCA,CBは全Pκで滑らかである。

Rのsmoothnessにも単に微分の下界を使わない。
全proper旗模型の二つの旗を同時に保持した閉incidenceを、R\Nの上へ制限する。
全66表によれば、この上の全旗対は最後のbig relative orbitに入る。
このorbitはlocally closedだから、全incidenceのこの上の部分において、
big-orbitのrank8線形交叉束は閉部分スキームとなる。
従ってそのR\Nへの射はproperである。

詳しくは、big orbitが閉となる旗積の開Uを選べる。
R\N上の全incidenceの幾何点は全てU内にあるので、そのUへの逆像は全体であり、
big orbitとのintersectionは閉となる。sourceの像はもともとRに入り、
固定したbig orbit上の交叉は一定rank8なので、このsourceは滑らかなP7-bundleの開部分である。

全Pκでその幾何fiberは一つで、20次minorはprojective sourceの全19接方向がinjectiveであることを示す。
同じ長さ1/Nakayama論証によって、sourceはPκ近傍の被約Rと同型となる。
よってRも全Pκで滑らかである。

## 3. κ=0でも二次obstructionは非零

既存の二つの実source曲線は全κで同じ一次方向vを持ち、その全二次差は−18 z²∂x∧∂yである。
κ≠0では既存の汎関数κe72−e63がtangentの和を零化し、差への値−18κが非零となる。

κ=0では別の汎関数

    L0=(5/22)e68+e72

が全tangentの和を零化し、二次差への値は−18となる。
68,72は同じcanonical150係数の0始まりindexである。
実際、κ=0でもtangentの和はrank37、二次差を加えるとrank38である。
選んだvもR接空間の外にある。
したがって第三次を追う前に、全cokernelを調べた二次obstructionだけで全κを覆える。

## 4. embedded構造を残さず通常二重環となる理由

CAは滑らかな28次元であり、Dのtangent dimensionは全κで20である。
CB方程式のrestrictionから独立な一次項を持つ八本を取り、形式的implicit function theoremで除くと、
Dは正則20次元環Sの商として書ける。
その被約台は滑らかなR19だから、S内でRのprime idealは一つのregular parameter(t)である。

Rに沿うgeneric primary lengthが2であることは、既存の全source二次obstructionの定理で既に証明した。
従ってDを定める残りのideal I⊂Sは、height1 prime(t)でlocalizeすると(t²)となる。
Sは正則なpower-series UFDなので、Iの全ての元はt²で割れ、I⊂(t²)である。
これにより、embedded構造を作るu*tのような一次のt項はgeneric length2と両立せず、最初から排除される。

各Pκでの非零二次obstructionは、Iのある元の二次項が非零な定数倍のt²となることを意味する。
その元はt²*unitなので、逆包含(t²)⊂Iも従う。
したがってI=(t²)であり、完成局所環は定理の通常二重環となる。
generic m=2だけでembedded primesがないと推論したのではなく、
全κでのregular source・smooth reduced support・非零quadratic unitを合わせている。

## 証明書と残る範囲

* `code/ambient_R19_open_double_20261103_v63.py`
* `reports/AMBIENT_R19_OPEN_DOUBLE_20261103_v63/report.json`

四つのconstant minorの全行・列、κ=0の全150-covector、obstruction、rankを保存した。
AF=0やt=0の境界の全scheme structure、および全域primary decompositionは本定理の対象外である。
