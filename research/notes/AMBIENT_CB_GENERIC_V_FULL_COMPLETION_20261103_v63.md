# CB generic V境界のactual全完成環

2026-11-03。全actual QQ(V) odd firstjet計算と全proper fiberの適用鎖に対する独立監査PASSを必須とするcompanionである。
Vを固定した一点の型からgeneric CBへ外挿しない。

## generic Vで全proper sourceを取る

P0での全CB proper flag fiberは一台長2であり、その周囲ではproper写像はfiniteで全fiber長2以下となる。
全polynomial critical V曲線はV=0でP0を通る。
従ってその逆像開集合はV曲線に非空開集合を定め、generic Vでは全proper fiber長が2以下である。
generic source Jacobianはrank28、sourceは滑らかな29次元なので、その全fiberは一台長2となる。
この時点ではfinite写像を得ただけで、global proper modelの未証明のbirationalityを使わない。

QQ(V)上の全150係数で、GL19方向とV曲線の微分はrank20であり、保存した二十次minorは

    −5103*V^15*(V²+1)/1562500

である。これらは全critical schemeに含まれる。
全source critical cotangentのrank9も確認したので、generic critical schemeは正則20次元となる。
GL orbitだけのrankは19であり、ここではmodulus Vの一方向を明示的に残している。

実28 target座標と一branch座標から、正則29次元のquadratic base Rを取る。
finite sourceとactual completedimageは

    B=R[t]/(t²−nu),       A=R⊕I*t ⊂ B

である。source critical quotientは B/(t,I)=R/(I,nu) なので、J=(I,nu)はRの正則九生成idealである。
全150 actual odd coefficientのQQ(V) firstjetはrank9である。
従ってIはJ/mJを全て張り、NakayamaによりI=Jとなる。
非零c∈Iに対しt=(ct)/c∈Frac(A)となるので、finite normal BはAと同じ分数体を持つ。
ここで初めてBをactual normalizationと同定する。

## 全完成型

Rのformal座標を選ぶとI=(nu,s1,...,s8)であり、generic Vでのfull cone completionは

    A=C[[g20,s8,t²,s1*t,...,s8*t,t³]]
      ⊂ B=C[[g20,s8,t]].

これはCAから型を移した推測ではなく、独立なCB QQ(V)計算から得られる。
formal coreが正則九生成odd idealであることを証明した後の純粋代数は、受理済みCA companionと同じである。
九odd target Y_i=b_i*t、b=(nu,s8)について、36 Koszul関係と45 pair-product関係が全kernelを与える。
pair productsでR+ΣR Y_iに還元し、残る関係を正則列bの第一Koszul syzygyで全て生成する。
この81関係は別のGBを仮定せずに完全であり、Macaulay2は本packetで実行しない。

coneではdim29、embedding dimension38、depth21、multiplicity2である。
GL E_yyのactual列は点P自身に等しく、scalar方向を除いたprojective版ではg19となり、dim28、embedding dimension37、depth20、multiplicity2となる。
normalizer conductorは(t²,s8)B、target conductorは(t²,s8,t³,s8*t)Aである。
conductor quotientは正則C[[g20]]、B/A=C[[g20]]*tであり、正確なformal singular、nonnormal、non-CM locusはこの滑らかなconductor supportである。
Aはdomainかつanalytically irreducibleである。
t²,t³∈Aによりnormalizationはelementary subintegralなfinite universal homeomorphismとなる。

これは**generic V**での結果である。
V=2など特定のsampleが全proper finite開集合に入るとは仮定せず、全Vの同型型を主張しない。
V=0、frameのV²=2除外、その他のspecial Vもこのgeneric定理で分類しない。
B18原点の全CA/CB odd ideal、三変数cubic、actual intersection Dの型も未解決のままである。
