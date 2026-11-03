# R19 に沿う実交叉の一般primary algebraは正確に二重点である

2026年11月3日。対象は二つの被約成分CA,CBの実scheme intersection D=CA∩CBである。
既にDred=N27∪R19、Rの一般点でCA,CBは滑らか、Dの接空間は20次元、Rは19次元と分かっている。
ここでは接空間の超過だけで長さを推定せず、全source表示の二次obstructionを計算する。

## 定理

K=C(R)をRの関数体とすると、DのRに沿うgeneric primary algebraは正確に

    K[epsilon]/(epsilon²)

である。従ってgeneric multiplicityは2であり、一般R成分はreducedではない。
これは全域primary decompositionや、κ=0を含む全特殊点での長さを決める主張ではない。

## 1. 一般横断の接埋込み次元は1

Rは一般点で滑らかであり、Nはその一般点を含まない。
従ってRのgeneric pointでDをlocalizeすると、residue field Kを持つArtinian local algebraとなる。
一般閉点のdim T(D)=20とdim R=19から、その横断cotangentの次元は1である。
標数零の等標数Artinian local algebraにはcoefficient field Kが取れ、極大idealが一つの元で生成されるので、
そのgeneric algebraはあるm≥2についてK[epsilon]/(epsilon^m)となる。
未決だったのはmである。

## 2. 実際の二つのsource曲線

一般正規形を

    Pκ=∂y∧(4x²∂x−2xz∂z+(2xw+h²+κhz+z²)∂w)

とする。
CAの固定旗sourceでは

    f=−4x²,
    H0=h²z+(κ/2)hz²+z³/3+2xzw

を使い、source変分deltaH=3z²wを取り、実際のGL5変換

    gA(epsilon)=I+2epsilon E_(x,z)

を適用する。E_(x,z)²=0なので逆行列はI−2epsilon E_(x,z)であり、打切り近似ではない。

CB側は、置換(1,3,2,0,4)を戻した元のsource座標でf0=6hxとする。
元の十九個のfiber係数を使い

    delta f=−21x²,
    delta H=xzw

を取り、source group変換は恒等とする。
CBの比率µの速度はこの選択では零である。
コードではg=µfをsourceで先に乗算し、二次のdeltaµ*deltaf項を含めてから零になることを確認した。
必要なmixed termを式から落としたものではない。

両曲線のcanonical150係数の一次微分は厳密に同じベクトルvとなる。
Rの八source方向と全GL5軌道方向はaffine coneでrank20、vを加えるとrank21なので、
vは被約Rの接方向ではない。

## 3. 全150係数の二次差とannihilator

全有限source式を展開すると、二曲線の差は正確に

    P_A(epsilon)−P_B(epsilon)
          =−18epsilon² z²∂x∧∂y

となる。三次以上のepsilon項も零であり、有限変換の全恒等式として確認した。
canonical150係数の0始まりindexでは、差の二次係数は72番だけが−18である。

CA,CBの全tangent像の和をWとする。
これは同じ150係数空間でQQ(κ)上rank37である。
次の明示線形汎関数はWを恒等的に零化する。

    Lκ=κ*(∂x∧∂y の z²係数)−(∂x∧∂y の hz係数).

二次差への値は正確に

    Lκ(P_A,2−P_B,2)=−18κ

である。κ=1で非零であり、一般κでも非零となる。
任意の二次source補正を加えても、二次像の変化はCAとCBの全一次tangentの和Wへ入るため、この値は変えられない。
したがってvを二次まで共通に持ち上げることはできない。

一般のgroup速度Mに対しても二次係数は

    Lie_M(deltaP)+(Lie_M²−Lie_(M²))(P0)/2+source_mixed

で計算しており、上の直接有限curve恒等式と全150係数で一致する。
quadratic Jacobi式の一部やnormalization tangentだけを比較したものではない。

## 4. m=2への帰結

CA,CBのproper sourceはこの一般R点の近傍で同型であり、そのsmooth source chartは全形式的変形を表す。
従って上のnonzero obstructionは実際のtarget intersectionに対するobstructionである。
また正規形族とPGL5作用はRを支配するので、κの非空な開でのこのnonzero classはRのgeneric pointで非零となる。

もしm≥3なら、K[epsilon]/(epsilon^m)からK[epsilon]/(epsilon³)への商により、
唯一の非零横断tangentは必ず二次まで持ち上がる。
しかし上の実source obstructionはそれを禁止する。
よってm≤2であり、既に非被約なのでm=2と決まる。

κ=0では今回の汎関数の値が零になるため、同じ計算だけではその特殊点の長さを決められない。
一般二重点という結論を、全Rが一様な二重構造だという主張へ広げない。

## 証明書

* `code/ambient_R19_generic_double_20261103_v63.py`
* `reports/AMBIENT_R19_GENERIC_DOUBLE_20261103_v63/report.json`

実sourceの二つの一次係数、全150二次差、全tangentを零化するcovector、mixed項、
直接有限curve恒等式、R接空間からの独立性をQQ(κ)上で保存・検査した。
