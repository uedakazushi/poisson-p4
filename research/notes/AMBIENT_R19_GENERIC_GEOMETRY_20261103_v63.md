# R19 の一般旗・有理性・実交叉の一般非被約性と全 CA/CB の R1

2026年11月3日。R19は全射影被約交叉N27∪R19の階数二成分である。
ここではRの一般点の幾何を扱い、t=0等の全特殊軌道は分類しない。
個々の成分のsmoothnessと、全Poisson空間のsmoothnessは区別する。

## 定理

1. Rは有理な19次元多様体であり、一般点のCA旗・CB旗は共に一意に戻る。
2. 全proper旗模型からCA,CBへの射はbirationalで、Rの一般点の近傍では同型である。
   特に個々のCA,CBはそこで滑らかな28次元多様体である。
3. 実交叉D=CA∩CBはRに沿って **generically nonreduced** である。
   一般閉点でdim T(CA)=dim T(CB)=28、dim T(D)=20だが、被約Rの次元は19である。
4. **全射影の被約CA,CBは共にSerreのR1条件を満たす。**
   これはrank4開だけの旧R1を、今回のrank2境界の議論で全域へ拡張する新しい結果である。

一般交叉の長さ・冪零次数・全primary structureはまだ計算していない。
全Poisson環の非CM性は別ノートのLPB包含によって扱う。

## 1. 正規形のintrinsicな三つの線

既存の一般正規形を

    Pκ=∂y∧Vκ,
    Vκ=4x²∂x−2xz∂z+(2xw+Q(h,z))∂w,
    Q=h²+κhz+z²

とする。canonical quotient fieldは

    C=(-8hx/5,12x²/5,-18xz/5,Q+2xw/5)

である。
全非零Qについてintrinsic factor L=C∂yは一意で、Cの全invariant linear hyperplanesは
h,x,zだけであり、そのラベルはcofactor比によってintrinsicに区別される。
これは `AMBIENT_R19_GENERIC_ORBITS_20261103_v63.md` の全割り算証明を使う。

またCの四成分は線形独立で、canonical二ベクトルには非零な定covector Casimirがない。
CAのAfとCBのBaの全点は二つの独立な定covector Casimirsを持つため、この点はいずれにも入らない。
従ってCAのproper持ち上げはAg上、CBの持ち上げはBb上に限られる。
どちらのsource constant factorも一意なLと一致する。

## 2. CAとCBの旗を区別するtrace identity

CAのAg表示をY=−f(ell,k)∂k+Hamiltonianと書く。
div Y=−f_kなので、canonical C=Y+(f_k/5)E_Uとなる。
marked line ellのcofactorをm=C(ell)/ellとすれば

    f=m*k−C(k),    f_k=5m.                  (CA)

CBのBb表示はradial項−f_k Ezw/2を持ち、div Y=−2f_kだから

    f=m*k−C(k),    f_k=(5/2)m.             (CB)

これらは商U=C5/L上の線形共役と倍率に不変な必要条件である。
ellはinvariant hyperplaneだからh,x,zのどれかであり、mは常に非零なxの倍数である。
m=f_kの倍数はbinary pencil Wに入るので、ell≠xならW=span(ell,x)でなければならない。

ell=h,k=xではf=−4x²、f_k=−8xで、CAの必要条件だけを満たす。
ell=z,k=xではf=−6x²、f_k=−12xで、CAの必要値−18xにもCBの必要値−9xにも合わない。
kをax+b ellへ取り替えてもこの不一致は変わらない。

ell=xの場合、CA条件は

    C(k)+(48/5)x*k∈C x²

である。この条件はQ≠0の下でk∈Cxしか許さず、独立性に反する。
CB条件は

    C(k)+(18/5)x*k∈C x²

である。k=a h+b x+c z+d wとすれば左辺は

    2a*hx+6b*x²+d*(Q+4xw)

なので、a=d=0で、W=span(x,z)が一意に定まる。
実際にこの旗でf=6xzのCB表示があり、全canonical150係数で一致する。

よって全幾何的source flagsは一意で

    CA: (L,W,H)=(C∂y,span(h,x),Ch),
    CB: (L,W,H)=(C∂y,span(x,z),Cx)

となる。CBの比率はb=0で、F≠0だからsource点も一意である。
固定旗の線形canonical写像の単射性から、fiber係数にも別の持ち上げはない。

## 3. 全source微分とpropernessからsmoothnessへ

元のcanonical150係数で次の全source微分を計算した。

* CA:固定旗の二十方向と全GL5軌道方向。
* CB:固定b=0の十九方向、µ方向(F=6hxに対するg=F)、全GL5軌道方向。

これは単に選んだsubfamilyの微分ではなく、滑らかなproper sourceの全接空間を張る。
affine coneでのsource次元は29であり、QQ(κ)上の実像微分も双方厳密にrank29である。
κ=1で非零な29次minorを保存した。projective化でsourceと像の接空間から同じscale一方向が除かれるので、
projective sourceの28方向にkernelはない。

前節により全幾何fiberは一つの点で、その点のfiber tangentが零なので、fiberは長さ1の被約一点である。
propernessとfiber dimension上半連続性により、この点の近傍で射は有限となる。
target局所環をA、有限source代数をBとすると、sourceはirreducibleで射はdominantだからA→Bは単射である。
B/m_A B=CなのでNakayamaによりBは1で生成され、A=Bとなる。
従ってsourceからtargetへの射はこの近傍で同型で、CA,CBは滑らかである。
この同型の非空開があるので、両proper射は全域でbirationalでもある。
source微分が大きいだけからtarget smoothnessを推測したものではない。

## 4. 実交叉はRに沿って一般に厚い

上の二つの全affine tangent像を同じ150座標内で並べると、QQ(κ)上のunion rankは **37** である。
双方29なのでaffine交叉tangentは29+29−37=21次元、projectiveでは20次元となる。
κ=0,1,2,3でも同じrankを正確に確認したが、一般結果にはQQ(κ)上の全行列を用いる。

前節のsmoothnessにより、これらは個々の成分の真の全tangent spacesであり、
normalizationからの小さい微分像と取り違えていない。
実scheme intersectionのtangentはこの二つの交叉に等しい。

全被約交叉はN∪Rで、RはNに含まれないので、Rの一般点の近傍の被約台はRだけである。
Rは一般点で滑らかな19次元多様体である。
もし実交叉がそこで被約なら、一般点ではtangent dimensionも19となるはずだが、実際は20である。
したがって交叉はRに沿ってgenerically nonreducedである。
ここから一般厚みの長さ2を推論してはいけない。現時点では正確なnilpotence orderを計算していない。

## 5. Rの有理性

全rank2交叉のtableでRを与えたpaired-flag orbitは12次元で、その上のlinear intersectionは8次元である。
従ってsourceはそのorbit上のP7-bundleで、19次元である。
前節の一般旗一意性から、Rへのdominant射はgenerically one-to-oneであり、標数零なのでbirationalである。

paired-flag orbit G/(P∩wPw^−1)は有理である。
実際G/P上のbig cellで射影は自明化し、そのfiberはflag variety内のparabolic orbitである。
後者は有限個のBruhat cellsで覆われる既約空間なので、一つの稠密なaffine cellを持つ。
P7-bundleはZariski局所自明だから、その全sourceも有理である。
従ってR19は有理である。

## 6. rank2境界の穴を埋めた全域R1

旧R1定理はrank4開でのみ成立していた。
今回の全proper模型によれば、CAのrank≤2境界はAg27とAf25、CBのそれはBb27とBa26である。
CA,CB自身は28次元だから、rank2境界のcodimension1部分はAgとBbだけである。

今回のR点はそれぞれAgとBbに入り、その点でCA,CBは滑らかである。
滑らかな開部分は各既約divisorを非空に切るので、AgとBbの一般点も滑らかである。
他のrank2境界はcodimension≥2であり、R1を妨げない。
rank4開の既存R1と合わせると、全射影の被約CA,CBがR1を満たす。
Nの全rank2divisorに対する同様の主張は本ノートの入力ではなく、別に検証する。

## 証明書

* `code/ambient_R19_generic_geometry_20261103_v63.py`
* `reports/AMBIENT_R19_GENERIC_GEOMETRY_20261103_v63/report.json`

全二つの150行source微分、QQ(κ)rank、κ=1minor、CB十九係数、必要旗trace equationsを保存した。
測定1.66秒。全Poisson空間の局所環、Rの全特殊点、交叉の一般長さは未解決である。
