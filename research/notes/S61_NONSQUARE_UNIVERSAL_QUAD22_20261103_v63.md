# 六quadricsの普遍22次元 quadratic-syzygy kernel

2026-11-03。元の六quadricsのquadratic cofactor map168×210を扱う。大きな行列や新pivotを展開しない。全体の上界146・generic146は既証である。ここでは常に存在する22本を実係数で構成し、known-unit開での独立性を直接証明する。これは周囲210方向における点別rank146下界とは異なる。

## 1. 固定sourceと新しい恒等式

元lightを `(x0,x1,x2,x3,x4,Z,W)=(v0,v1,v2,v4,v7,3v0-v1+v6,v3)` とする。hZ=ell*Z,hW=ell*Wと書く。六quadのbasisを `(q0,q1,q3,q5,hZ,hW)` へ変える行列のdeterminantは **−80**、係数の分母はaだけである。この実変換を元sourceから照合した。

Z=W=0上のLは `L_x1=-40C`、

`q1bar=L²/8100`, `q3bar=L*(eL-400a*x0)/8100`。

LをZ,Wを含まない線形式としてそのまま持ち上げる。cubic

`(eL-400a*x0)*q1-L*q3`

は(Z,W)に入る。各単項式でZがあればZを先に一つ取り、なければWを取るという固定規則でactual quadratics Tz,Twを作り、

`(eL-400a*x0)*q1-L*q3=Z*Tz+W*Tw`

をliteral QQ polynomial identityとして検査した。従って六quadに対する新しいquadratic-cofactor syzygyは

`(0, ell*(eL-400a*x0), -ell*L, 0, -Tz, -Tw)`。

その全係数をcertificate.jsonに保存する。

## 2. 古い21方向

六quad間の15Koszulに加え、universal linear relation

`W*hZ-Z*hW=0`

をlightで掛ける。ellのx1係数は `-16*C*J/243` という既知unitなので、x1を除いた六light倍を選ぶ。

六quadricsは全actual開で定数上独立である。実際、Z=W=0上で四qの定数関係があれば、さらにL=x0=0へ制限し、既証のq0/q5 pair rank2より両係数が零。残りはLとx0の独立性によりq1/q3係数も零となる。hZ,hWを含む六式の関係も、同じ制限の後はell*(bZ+cW)=0から零である。

よってKoszulとlinear-relation倍との重なりは正確にell倍のrelationだけ。具体的に、Koszulの定数結合とlinear-relation*Vが零なら、四q成分と六quadの定数独立性から、KoszulはhZ,hW間の一つに限られ、Vはellの定数倍となる。選んだ六lightのspanはellを含まないので、上の21本は独立である。

## 3. 新しい22本目の独立性

Z=W=0上で

`Delta=ellbar-(2J/1215)*L`

と置く。actual coefficientsはDelta_x1=0および

`3*Delta_x2-(7/2)*Delta_x4=-80*a*d²/81`

を満たす。右辺はad!=0上のunit。L_x1=-40Cなので、ellbarとLは全adCJ開で比例しない。

新syzygyが古い21方向に属すると仮定する。Z=W=0へ制限する。q0/q5成分が零であることと四qの定数独立性から、残るKoszulはq1,q3間のものだけとなる。新しいq3成分は−ellbar*L、Koszul側は定数kに対して−k*L²/8100。Lを消去して `ellbar=k*L/8100` を得るが、上のknown-unit identityに反する。従って新しい方向も独立である。

以上、全 `a*d*C*J!=0`、特にFULLUでは22本のquadratic-cofactor kernel directionsが点別に独立である。新しいpivot探索、未知因子の除去、有限体は使わない。

## 4. 加群kernelと周囲rankの区別

22本はR168内のsubbundleをなし、quotientはlocally free rank146である。既証generic rank146により、quotientからR210への誘導写像はgenericに単射。baseはintegral、sourceはlocally freeなので、その加群kernelは零となる。従ってこの22本は**加群としてのkernelの全basis**であり、quadratic-product row moduleはlocally free rank146となる。

この最後の写像は特殊化後に単射とは限らない。従ってambient cokernelのflatness、全点rank146、fullquartic全射は結論しない。生成source側の可逆縮約と、周囲への点別injectivityを混同しない。global free146 basisや具体的146pivotも主張しない。

再現は `python code/s61_nonsquare_universal_quad22_20261103_v63.py --output <new-output>`。actual sixquad、basis変換、liftした全syzygy、二つのknown-unit anchorsを約1秒で照合し、cofactor polynomial matrix168列は展開しない。
