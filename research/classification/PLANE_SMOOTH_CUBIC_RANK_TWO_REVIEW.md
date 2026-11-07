# smooth cubic零点を持つrank二plane：全normal-weight源の独立検証

2026年11月7日。主義務 `E_reduced_nonnormal`（一般底Gの残余排除）。
**判定：記載された仮定でPASS。**
S=P²がPoisson不変でS⊂SingD、π|Sが一般rank二、その零divisorがsmooth irreducible cubicなら、
rank四πの実退化因子Dは既約・被約quinticにならない。
横断ordinary node、normalization/conductorのsmoothness、vertical divisor不存在は追加しない。

独立実装は `code/classification_plane_smooth_cubic_rank_two_review.py`、
証明書は `reports/CLASSIFICATION_PLANE_SMOOTH_CUBIC_RANK_TWO_REVIEW_20261107.json`、
有限契約は `jobs/plane_smooth_cubic_rank_two_review.json`。
producer・CAS・repository arithmeticをimportせず、standard-library Fractionと疎多項式のみを用いた。
全ternary cubicの十係数、全weight一の三十係数、全weight二の九係数をsymbolicのまま保持する。
特定のHesse cubic、数値a、有限個のsmooth parameter標本から一般化していない。
最終契約gate PASSの後の実行は約0.266秒、最大RSS13,988 KiB（CPU170秒、memory1,536 MiB上限）。
gate自体は数学の証明ではない。

## 1. 接続の既済入力とscalar部分を固定する橋

`PLANE_SMOOTH_CUBIC_CONORMAL_CONNECTION_REVIEW.md` の別作者による幾何監査を入力とする。
そこで受理したのはconormal接続のtracefree部分零であり、元のO(−1)のscalar接続零ではない。
この独立reviewは、そのmonodromy証明を重複して再計算しない。

base座標x=(x₀,x₁,x₂)、normal座標z,wとし、smooth cubicをfと書く。
canonical lift

\[
 P_f=f_0\partial_1\wedge\partial_2+
 f_1\partial_2\wedge\partial_0+f_2\partial_0\wedge\partial_1
\]

の誘導するconormal接続を基準に取る。P_fはstrict Poissonで、同じplane sectionを与え、接続はscalarである。
元の接続との差はα·id、α∈H⁰(P²,Ω¹(R))となる。
両接続のflatnessからdα=0である。
局所R=(t=0)で一位極を持つ閉一形式を書くと、接方向のdt∧dy/t²係数から
接方向poleの係数はtで割れ、αはlogarithmicになる。
既済のresidue完全列によるH⁰(Ω¹(logR))=0からα=0である。
従って全conormal接続はcanonical接続に一致するが、scalar接続自体を零としたのではない。

## 2. 全150係数からnormal filtrationと完全商を再構成

normal dilationで係数のnormal次数からnormal vector数を引いたweightを使う。
全二次斉次bivector150係数のraw weight寸法は

| weight | −2 | −1 | 0 | 1 | 2 |
|---|---:|---:|---:|---:|---:|
| raw係数 | 6 | 42 | 57 | 36 | 9 |

である。全25個のEuler外積のrank24も独立に生成した。
x₀≠0で各ratio座標の括弧を

\[
 \overline\Pi^{ij}=\Pi^{ij}-x_j\Pi^{i0}-x_i\Pi^{0j}\quad(x_0=1)
\]

と直接計算した。mixed括弧のnormal次数零、normal-normal括弧の次数零と一を全て消す線形条件は、
150列に対してrank42、kernel108である。
そのkernelは非負raw102とweight−1のEuler六方向の和で完全に張られる。
従ってEuler商後に負weightは存在しない。
これはinvarianceと{I,I}⊂I²を課した完全projective sourceであり、chosen ansatzではない。

weight零raw57上で、base restrictionとmixed normal一次接続が共に零となる全条件を
同じratio括弧から作るとrank42、kernel15を得る。
kernelはweight零Euler十二方向とnormal-normalのbinary quadratic三方向の直和で尽くされる。
全chart等式はdense affine chart上の多項式恒等式なので大域sectionの零と同値であり、
未処理のx₀=0 chartを捨てていない。
§1のcanonical接続との一致により、全weight零代表は正確に

\[
 P_0=P_f+K(z,w)\partial_z\wedge\partial_w
\]

へ取れる。

weight一raw36からEuler∧(zV+wW)の六方向を商すると、次の三十係数源が一意に残る。

\[
 P_1=zL+wM+z^2D\wedge\partial_w
       +zw(B\wedge\partial_z+E\wedge\partial_w)+w^2C\wedge\partial_z.
\]

L,Mはbaseのlinear coefficient bivectorで各九係数、B,C,D,Eはconstant base vectorで各三係数。
この三十列のrank30、Euler六列のrank6、結合rank36を全150座標で確認した。
weight二にはEuler方向が無く、完全源は

\[
 P_2=z^2A+zwA'+w^2A''
\]

の九係数（A,A′,A″はconstant base bivector）だけである。

## 3. 全smooth fに有効な唯一のinjectivity

constant base vectorXについてU_X=[P_f,X]=−P_(Xf)と書く。
Xf=0なら、X≠0の場合は線形座標変更でfが二変数だけの三次式になる。
そのplane cubicは残りのcoordinate pointで特異となるのでsmoothnessに反する。
従ってXf=0⇒X=0、またU_X=0⇒X=0である。
全cubic十係数のまま、Euler_base∧U_X=−2Xfを独立に検証した。

linear base bivectorLを三次元volumeでlinear vector lと同一視すると、
[P_f,L]の唯一のstrict三vector係数は符号を除きgrad(f)·curl(l)である。
この全九係数の恒等式も直接確認した。
curl(l)はconstantだから、closedならcurl(l)=0。
curl mapのrank3、kernel6も確認したので、Lはquadratic gのJacobian bivectorP_gである。

## 4. K=azw、a≠0：共鳴の無い完全kernel六

全bilinear SchoutenからEuler外積を取り、五つの四添字・全70単項式による全350係数を構成した。
二base方向とz,wを持つ成分のz,w次数三の係数は、base bivectorとして正確に次となる。

| normal monomial | coefficient |
|---|---|
| z³ | U_D+a Euler_base∧D |
| z²w | U_E−aL |
| zw² | −U_B−aM |
| w³ | −U_C+a Euler_base∧C |

さらに三base方向とwを持つ成分のz²係数はEuler_base∧U_D、
三base方向とzを持つ成分のw²係数はEuler_base∧U_Cである。
§3からD=C=0となり、四係数表は

\[
 L=a^{-1}U_E,\qquad M=-a^{-1}U_B
\]

を強制する。これは全a≠0・全smooth fで有効であり、固有値の共鳴値を除外しない。
E=aY、B=−aZと置くとP₁は実normal→base shear X=zY+wZの軌道微分そのものである。
この全六係数解を残る全Jacobi成分へ戻し、全350式が恒等零となることも確認した。

## 5. K=az²、a≠0：重根型の完全kernel六

同じ全350式の二base・二normal成分は正確に

| normal monomial | coefficient |
|---|---|
| z³ | U_D−aL−a Euler_base∧E |
| z²w | U_E−aM+a Euler_base∧B |
| zw² | −U_B+2a Euler_base∧C |
| w³ | −U_C |

である。最後からC=0、B=0、M=a⁻¹U_E、L=a⁻¹U_D−Euler_base∧Eとなる。
D=aY、E=aZと置けば、shear X=zY+wZの軌道微分から
Euler∧(azZ)だけを引いたcanonical代表と正確に一致する。
全三十係数kernelがこの六係数で尽くされること、および全350への再代入を独立に確認した。
二つのnonzero Kは複素数体上のbinary quadraticの全GL₂軌道を尽くす。

## 6. 真の線形shearとweight二の全injectivity

X=zY+wZのlinear matrixの平方は零なのでg=Id+Xは可逆で、inverseはId−Xである。
独立算術では、x→x−zY−wZの全係数代入と全vector方向のpushforwardを行い、
g_*P₀のweight零がP₀、weight一が上記実orbit微分、残りがweight二以下であることを全変数で検証した。
単なるtangent vectorを実変換と同一視していない。

元のπ=P₀+P₁+P₂にinverse shearを施すと、weight一が零となり、
global quadratic degreeから残りはP₀+P₂′だけである。
[P_f,P₂′]=0も全九係数・全f十係数で確認した。
weight二Jacobiの二base・二normal成分は、各base bivector係数Qについて正確に−2KQとなる。
K≠0は多項式環の零因子でないためP₂′=0となる。
従ってπは実PGL変換でP_f+Kへ戻る。
符号付きEuler Pfaffianは全f十係数のまま

\[
 D_{P_f+K}=3f(x)K(z,w)
\]

と一致する。相異二根は可約、重根は非被約なので、既約被約quinticを与えない。

## 7. K=0の全補集合と保持する射程

K=0では二base・二normal成分は
z³U_D+z²wU_E−zw²U_B−w³U_Cであり、§3からB=C=D=E=0。
残るlinear JacobiはL=P_g、M=P_h（g,h quadratic）である。
このkernelはweight一の十二次元であり、P₂は九係数全部がbase-base方向のみ。
よって全π=P_f+zP_g+wP_h+P₂は三つのbase vector方向に含まれ、射影rankは二以下。
higher-weight Jacobiを全て解く必要なくrank四不存在が従う。
K=0をnonzero Kの正規化やa反転に吸収していない。

これでsmooth irreducible cubic零divisorを持つrank二planeの全源が排除された。
singular/reducible/nonreducedな零cubic、rank零plane、plane以外の曲面は本reviewの補集合である。
既約非正規quintic全体又はG全体の尽くし、新componentのE帰属は主張しない。
