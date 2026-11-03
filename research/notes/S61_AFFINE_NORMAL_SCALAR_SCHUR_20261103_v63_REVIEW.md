# Scalar Fittingと15 Schur chartの独立監査

2026-11-03。`S61_AFFINE_NORMAL_SCALAR_FITTING` と `S61_AFFINE_NORMAL_SCHUR_COVER` を承認する。
全45一次式が**scalar entry**のK定数spanに入ることと、全未知Bを覆う15個のscalar idealへの正確な移行を区別する。
Aのparameter annihilation、具体的Lambda、chart unit、尽くしは主張しない。

## 1. 全3350 QQ entryと構造上界

全134 remainderから別のcompanion actionで全670 vectorを構成した。
さらに元S61の全134式のq二次項を実T01 betaへ直接代入し、保存E frameでL_i(n,S1,S2)へ移した。
centeringで二次項が変わらないことを使い、作者のload_fiber/Q_i構成は呼ばない。
全3350 entryで、Hを含まない二次部分が厳密にdelta_(j,k) Q_i(L)であることをQQ全係数で確認した。
LはHを含まず、他の二次部分は全てH*x_aに入る。

generic betaでも、元式はqに二次、q=P(h)+L、Pのh次数≤4である。
non-H二次Q_i(L)はh^0にだけ現れ、h^j倍(j≤4)で割り返されない。
残るdivision由来の二次はH*xだけなので、二次spanの上界134+45=179は普遍的である。
凍結点のmodp rankや欠落六monomialから一般上界を外挿していない。

全1035二次と45一次の1080列を保持し、全237811 native entry、QQ分母unit、選択224行列のmonomial labelを照合した。
別C++の整数modular Gaussian算法で、選択224次determinantは21497 mod32003と一致した。
generic係数は親quinticの承認済みbeta有理recipeでspecializeできるため、この非零minorはQQ/generic rank下界224を与える。
構造上界179+45=224と合わせて全rank224、二次rank179となる。
二次射影のkernelは45次元で、全一次空間に等しい。224minor逆行列による有限cofactor recipeも正しい。
巨大QQ逆行列や45本の展開cofactorを計算済みとは主張しない。

## 2. 実moduleとFittingの意味

formal R'=K[n42,S1,S2,H]上のpresentationは実Aを保つ。
実際、Hをg(h)で消すと、(h,S1,S2,n)から元qへのmapは
z1=p1(h)+S1,z2=p2(h)+S2と可逆E frameの合成であり、三角的な多項式座標同型である。
元134式が全て残るので、自由placeholderによる余計なschemeを導入していない。
normal polynomial recipesはこの商で成立する恒等式である。

scalar entriesは全てm=(n42,S1,S2,H)に入り、その定数spanが全linearsを含むのでFitt4(A)=mである。
これらは異なるtarget位置のentriesの結合であり、vector module rowの結合ではない。
mがAをannihilateすることも、非線形42-normal basechange後の同じ定数spanも結論しない。

既承認のknown-five split・normal閉埋込みにより、R'→Bは全射で、
Aはmoduleとして(R'/m)^5 directsum R'/I_Bである。E5の代数がK^5にsplitするとは仮定しない。
Bはmを避けるためm+I_B=R'であり、各点での直接和Fitting計算から

    Fitt_j(A)=m^(5-j),  j=1,...,4,
    Fitt0(A)=m^5 intersect I_B

を得る。これはI_Bのnilpotentsを落としたradical公式ではない。

## 3. 非被約B上でも固定bottom四行が使える

R'→Bが全射でknown-fiveとsupportが交わらないので、presentationのright-exact basechangeは
coker(Phi tensor B)=Bを与える。quotient rowはtau=(1,h,h²,h³,h⁴)である。
tauのkernelは列(e_j-h^j e0),j=1,...,4で張られ、bottom四座標への射影が逆となるfree B4である。
この恒等式を記号hのまま確認した。AのflatnessやBの被約性を使わない。
従ってPhiのbottom四行はB上でfree B4へ全射である。

670×6のLambdaを4020個の独立constantとしL=K(lambda)とする。
Bの各幾何点でLambdaから4×6行列へのmapは全射線形であり、rank≤3条件は余次元3である。
dim B≤2は既証の実finite two-coordinate projectionから来る。
bad incidenceの次元は4020-24+21+2=4019以下で、parameter次元4020より小さい。
従ってgeneric constant extension上ではrank落ちはBの全点で起こらず、15個のbottom4×4 minorsが全B_Lを覆う。
非proper射影でも像の閉包の次元上界で十分であり、generic点だけを検査した主張ではない。
nilpotentsは開被覆のexact localizationと有限moduleのNakayamaによって保持される。

## 4. Schur scalar idealは実schemeを保つ

deltaを一つのbottom4 minorとする。originでは全Phi entryが零なのでdeltaも零で、known-fiveはこのchartに存在しない。
四本のpivot columnは元670列のconstant結合であり、追加してもimageは変わらない。
deltaを逆にして四targetを消すと、残るtargetは元の1である。
各旧columnの関係はg_j/deltaで、

    g_j=det[pivot4columns,oldcolumn_j]

となる。普遍5×5 determinantとadjugateの式を全記号係数で確認した。
従ってcokernelはmoduleとしてだけでなく、1をgeneratorとする実algebraとして
R'_L[delta^-1]/(g_0,...,g_669)に等しい。nilpotentsを含めて正確である。
次数はdelta≤8、g_j≤10である。

15個全てのlocalized idealがunitであることはB_L=0と同値で、L/Kの忠実平坦性によりB=0とも同値になる。
具体的なQQ Lambdaを選んだ証明書、いずれかのunit、全未知の空性は計算していない。
新しい次数probe・Macaulay・GBは行っていない。

review code: `code/s61_affine_scalar_schur_review_20261103_v63.{py,cpp}`。
receipt: `reports/S61_AFFINE_NORMAL_SCALAR_SCHUR_20261103_v63_REVIEW/`。
