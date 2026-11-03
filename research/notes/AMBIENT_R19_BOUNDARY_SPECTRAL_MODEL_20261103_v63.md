# R19 の proper 二次スペクトル模型と、全被約境界 R∩N=B18

2026年11月3日。ここでN,Rは全射影CA/CB交叉の既存の被約成分N27,R19である。
本稿の交叉等式は全て **underlying reduced closed varieties** の等式である。
R∩Nのscheme厚みや全Poisson空間のprimary decompositionは決めない。

## 定理

Rには、滑らかな9次元旗多様体上の相対quadricからなる **normalなprojective19次元proper birational模型** がある。
その一つのirreducible Cartier divisorのproper像がRの全被約境界であり、

    (R∩N)_red=B18

は既約18次元である。

B18は、座標(h,x,y,z,w)における九係数の全線形族

    P=∂y∧D,
    D(h)=a*x², D(x)=0,
    D(z)=b*hx+c*x²,
    D(w)=G2(h,x,z)

のPGL5飽和閉包である。G2は任意の二次式で六係数を持つ。
この族からのfull-flag P8-bundleもproperにB18を覆う。

さらにB18の中に閉既約部分B17,B16があり、

    B16 ⊊ B17 ⊊ B18

となる。次元は添字の通りである。
元のbig R8 cellでt=q4=0とした全像はB16に等しいので、その像だけではR∩Nを尽くさない。
Nの二つのrank2境界をNg26,Nf25と書けば

    (R∩Ng)_red=B18,
    (R∩Nf)_red=B17.

また被約Rは全域でR1を満たす。ただしproper模型のnormalityから全Rのnormalityを結論しない。

## 1. 部分旗と十二次元canonical線形空間

底は

    L⊂C5:直線,
    X⊂K⊂ann(L): dim X=1, dim K=3

という旗である。双対ではF(1,2,4;5)で、次元は9である。
基準旗をL=C∂y,X=Cx,K=span(h,x,z)とする。
canonical bivectorをP=∂y∧Dとし、

    D(h)=x*((r−13t)h+b*x+c*z),
    D(x)=12t*x²,
    D(z)=x*(d*h+e*x+(−r−13t)z),
    D(w)=2t*xw+G2(h,x,z)

と置く。十二変数は(t,r,b,c,d,e,G2の六係数)である。
全てdivergence zeroであり、Pは既にcanonical代表である。
P=∂y∧DかつDはyに依存しないので、全十二次元線形空間が生のPoisson条件も満たす。

この部分旗のGL5 parabolicは

    row x: xのみ,
    row h,z: h,x,zのみ,
    column y: row yのみ

を保つ16次元連結群である。
全16 Lie基底と十二個のcanonical基底について、192列の全150係数を計算し、同じ十二空間へ戻ることを確認した。
従って全群で不変であり、rank12の随伴ベクトル束が定まる。
canonical写像のrank12も検査済みで、projective mapにbasepointはない。

## 2. 一つの相対二次式がproper R模型を与える

次を課す。

    q=r²+c*d−25t²=0.

parabolic Lie作用において、tのcharacterはM_yy−M_xx、qのcharacterはその二倍である。
全係数で確認したため、これはprojective bundle上のwell-definedな相対二次式となる。
qのrankは4であり、fiberはP11内のirreducible quadric hypersurface、次元10である。
底の9方向と合わせ、sourceの次元は19となる。
sourceはhypersurfaceでS2、特異集合はfiberのvertex P7、従って全source内でcodimension3なのでR1でもある。
よってsourceはnormalである。

K上のoperator TをD(ell)=x*T(ell)で定める。
その行列は(h,x,z)で

    [[r−13t,b,c],[0,12t,0],[d,e,−r−13t]]

であり、q=0の下でcharacteristic polynomialは厳密に

    (lambda+8t)(lambda−12t)(lambda+18t)

となる。
t≠0では三固有値は相異なり、Xを保つKのbasisでh,zもeigenvectorsにできる。
canonical fieldからV=(D+8t*x*E_U)/5を戻すと

    V=4t*x²∂x−2t*xz∂z+(2t*xw+G2/5)∂w

となる。既存の有限w-shiftsはhx,x²,xz項を全て除き、Rの標準形V_Qを与える。
Qは零でもよく、これら全標準形がR\Nに入ることは既に証明済みである。

逆に全R\Nの標準形はt=1,r=5,b=c=d=e=0,G2=5Qとしてsourceに入る。
t≠0はsourceで稠密、sourceはprojectiveなので全閉像は正確にRである。
L,X,Kはこの開でintrinsicに戻るため、射はbirationalである。

## 3. t=0は全てを覆う一つのirreducible境界

sourceでt=0とすると、二次式はr²+c*d=0となる。
これはrank3のirreducible quadricで、その相対projective hypersurfaceの全次元は18である。
t=0はsourceのirreducible Cartier divisorである。

このときTはnilpotentでX⊂ker Tである。
三次元nilpotent linear algebraにより、必ず二平面Wを選び

    X⊂W⊂K,  T(W)⊂X,  T(K)⊂W

とできる。
regular nilpotentならW=ker T²、rank1の場合もimage TとXを含む二平面を選べばよい。
T=0なら任意である。
適合するbasis(x,h,z)を取ると、Dは定理の九係数chain族となる。
逆に全chain族はこのnilpotent条件を満たすので、t=0の全像はchain族の全proper像に等しい。

chain族はfull flag(L,X,W,K)の10次元底上のrank9線形空間である。
GL5 Borelの全15 Lie基底、135列でこの九空間の不変性を確認した。
従ってP8-bundleからの全射proper射がある。
指定点a=b=c=1、G2係数1,…,6で実affine像微分は厳密にrank19であり、
source上界18と合わせて像Bの次元は18となる。

chain族は全てNgに入る。
実際、N側のmarked coordinateをh_N=x、もう一つをk_N=hとすれば

    f_N=−a*h_N²,
    H_N=−(b*h_N*k_N+c*h_N²)w+integral G2(k_N,h_N,z) dz

というCAのN_g表示に一致する。全行列でこの恒等式を検査した。

任意のP∈R∩Nはproper sourceへ持ち上がる。
t≠0のsource像は全てR\Nだから、その持ち上げはt=0に限られる。
逆にt=0の全像はNg⊂Nに入る。
従って **全被約R∩Nが正確にB18** であり、他の旗退化の像は残らない。

## 4. 十七次元と十六次元の内在的な部分

B17を次のtriangular Hamiltonian族の全proper像とする。

    D(h)=D(x)=0,
    D(z)=−f2(h,x),
    D(w)=G2(h,x,z).

底は(L,W2,K3)、dim=9であり、rank9線形束のP8からのproper像である。
全144 Lie列、binary f2 discriminantのrelative characterを検査した。
実affine像rank18により次元17である。
B16はf2が平方となる閉部分のproper像で、irreducible source上界16、affine像rank17から次元16である。

chain族でa=0とするとf2=x*(-b*h−c*x)であり、C上の任意のbinary quadraticはlinear factorを持つため、
全G飽和後の像は正確にB17となる。
b=0ではD(h),D(z)は共通のx²を持つ同じ定方向であり、h,zを線形変換すれば平方f2型へ移る。
よってその全像はB16である。
元のbig R8でt=0の式はD(z)=−q0x²、G2任意であり、全G飽和は正確にB16である。
これで単なるt=0式だけでは全B18を見落とすことが分かる。

B17⊂NfはHamiltonian H=f2*w+integral G2 dzで明らかである。
逆にBのchain表示でNfへ入る点は少なくとも二つの定covector Casimirsを持つ。
a=0またはb=0なら上のB17へ入る。
a,b≠0なら、このCasimir条件はG2∈span(x²,hx)を強制する。
G2=mu*hx+nu*x²の場合

    w'=w−(mu/b)z−((nu−mu*c/b)/a)h

はCasimirとなり、残るfieldは平方f2型でB16に入る。
したがって(B∩Nf)_red=B17であり、(R∩Nf)_red=B17となる。
B⊂Ngなので(R∩Ng)_red=Bである。

## 5. Rの全域R1

chainの指定点a=b=c=1,G2係数1,…,6では、Lは一意なintrinsic factorである。
xは唯一の定covector Casimirで、Kは「D(ell)がxで割れる線形形式」の全空間として一意に戻る。
G2 mod xが非零なのでK=span(h,x,z)である。
t≠0のsourceはNへ入らないため、この点の全source fiberはt=0上のこの一つの旗に尽きる。

source quadricのgradientはc方向に1で、sourceはこの点で滑らかである。
quadric tangentの全11方向とGL5作用のcanonical微分はaffine rank20を持ち、projective sourceの全19接方向にkernelはない。
従ってfiberは長さ1であり、proper有限化・Nakayamaによりsourceはtarget Rの近傍と同型となる。
この点でRは滑らかである。

R\Nは既に全域で滑らかであり、残る境界は既約divisor Bだけである。
Bを切る非空なsmooth openがあるのでBのgeneric pointも正則であり、全RはR1となる。
sourceがnormalという事実だけからRのnormalityを結論したのではない。

## 証明書と残るscheme問題

* `code/ambient_R19_boundary_spectral_model_20261103_v63.py`
* `reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63/report.json`

全192+135+144 Lie列、quadricとtのcharacters、characteristic polynomial、
B18/B17/B16の全image minors、Rのboundary immersion minor、正確なN表示をQQで確認した。
測定7.73秒。report内ではB17,B16をT17,T16とも表記している。
R∩Nのscheme厚み、B自身の全局所環、全Poissonのglobal primary decompositionは未解決である。
