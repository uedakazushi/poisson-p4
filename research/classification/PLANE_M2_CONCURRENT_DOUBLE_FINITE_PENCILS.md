# concurrent三線・m2 double：全残源は十七のwhole a-linesへ

2026-11-08。主義務G/E。実D既約被約quintic、rank二Poisson plane、
reduced concurrent zero cubic f=xy(x-y)、normal位数二、generic double normal quadratic。
全moving rootは既空。constant rootのh=f、L^2M、L^3を扱う。
本稿はpureproofで、新しいmatrix/profile/GB/jetを実行しない。
十七族の全higher支持判定は次の別義務であり、ここでは全G完了を主張しない。

## 1. 原全normalformsとratio-Hamiltonian

`PLANE_M2_DOUBLE_ROOT_CONSTANT_NORMAL_FORMS.md` 第4–7節の全sourceを使う。
normal rootu、Puv=k u^2、X-Y=Delta、tracefree X,Y,Zとする。
全raw36P1/全九P2でRu=0、Rv=-2k u。
canonical strict liftのHamiltonian M=Pi sharp d(v/u)は

\[
 M(u,v)=\tfrac{3k}{5}(u,v),\qquad
 M(b)=(A(q)-\tfrac{2k}{5}I)b+B(u,v),\quad q=v/u,
\]

\[
 A(q)=Z-q\Delta.\tag{1}
\]

Bは全P1 mixedから来るbase-independent rational homogeneous degree一vector。
全higher係数を零にせず、全P2も許す。
qはM-invariant。M(F)=0はcanonical strict Casimirから従う。
式(1)とvertical補正-2kE/5の全証明は
`PLANE_M2_DOUBLE_EQUAL_DIAGONAL_NIL_OFF_EXCLUSION.md` 第2節と同じで、
本稿ではDelta非零も保持する。

## 2. affine translationとnormal gradingの一般判定

K=C(q)上でkI-A(q)がinvertibleなら
T=-(kI-A(q))^-1 B、b'=b+Tで全translationを消せる。
逆はb=b'-T、normal total degree一であり、source parameter以外のunknownを反転しない。
これは有理恒等式を読むK-affine changeであり、global PGLと呼ばない。

全体degree五Fのnormal total degree m summandは、base degree d=5-m、0<=m<=5であり、
Mの対応operatorは

\[
 A(q)\vert_{\operatorname{Sym}^{5-m}}+k(m-2)I.\tag{2}
\]

原Fはnormal位数二なのでm0/1は零。
leadingF2=u^2 hでA(q)h=0である。
m3/4/5の全operatorがinvertibleならFはexactly u^2 h(b+T)。
hがL^2M又はL^3ならK[x,y,z]上で二つ以上のpositive-degree factorsを持つため、
polynomial actualFはGauss lemmaにより既約になれない。
従ってoperatorのzero spectral wallsだけが指定sourceの必要残余となる。

## 3. h=fは全空

h=fではk nonzero、Delta=0、Z=(cx+dy)partial_zはnilpotent。
第2節の全operatorはm not equal twoでinvertible。
既 `PLANE_M2_DOUBLE_EQUAL_DIAGONAL_NIL_OFF_EXCLUSION.md` の全raw36/9証明により、
F=u^2f(b+T)はfraction field上三因子へ分解する。
よって全h=f sourceは空。Yの全二period及びZの全二periodを保つ結論である。

## 4. h=x^2y：全五b壁

全sourceは

\[
 Y=(ax+by)\partial_z,\quad X=Y+3y\partial_z,
 \quad Z=\operatorname{diag}(-1,2,-1)+c x\partial_z,
 \quad k=b+4,
\]

domain a+b+4 nonzeroである。
actualG2は3(a+b+4)x^2yu^2。
A(q)はtriangularでdiag(-1,2,-1)、off-diagonal c x partial_z-3q y partial_z。
従ってtranslation determinantは(k+1)^2(k-2)。
degree dのbase spectrumは3j-d、0<=j<=d。
第2節のhigher operatorsが零eigenを持つ全条件は

| normal m | base d | 全k |
|---|---:|---|
| 3 | 2 | 2,-1,-4 |
| 4 | 1 | 1/2,-1 |
| 5 | 0 | 0 |

であり、translation failureも既にこの表へ含まれる。
ゆえに指定sourceは全て

\[
 k\in\{2,-1,-4,1/2,0\},\quad
 b\in\{-2,-5,-8,-7/2,-4\}
\]

へ還元される。k=0はunitと仮定せず、そのwhole wallを残す。

更にcは0又は1へ正規化出来る。
weights(x,y,z,u,v)=(0,0,1,0,1)のactual GL5で、P0の全項はcharacter-1、
c x partial_zのoff混合項だけcharacter-2。
全体projective scaleを合わせるとa,b,k及びdiagonal Zは固定され、cだけがnonzero characterでscaleする。
c nonzeroならその値を一へ、c zeroなら閉補集合として零へ保つ。
inverseは逆nonzero diagonal scaleであり、higher coefficientsも全て同時に戻す。

これで全h=x^2y sourceは**五b壁 x 二c型の十whole a-lines**へ入る。
domainはa+k nonzeroで、そこ以外のaをgenericと仮定しない。

## 5. h=x^3：全七b壁

全sourceは

\[
 Y=(ax+by)\partial_z,\quad
 X=((a-3)x+(b+6)y)\partial_z,
\]

\[
 Z=(x+3E y)\partial_y-3E z\partial_z,
 \quad E=-\frac{b+8}{3(a+b+4)},\quad k=E(b+8).
\]

domainはa not equal four、a+b+4 nonzero。
これらは元actualG2=(a-4)x^3u^2と全normalizer equationsからの条件であり、
欠けた有効chartを分母で捨てたものではない。
b=-8ではE=k=0、Z=x partial_yというwhole boundary sourceを保持する。

Delta=(-3x+6y)partial_zなのでA(q)はtriangularでdiag(0,3E,-3E)。
E nonzero上、degree二spectrumは3E*j、j=-2,-1,0,1,2。
degree一spectrumは3E*j、j=-1,0,1。
operator(2)の全higher wallsは

| normal m | 全(b+8) |
|---|---|
| 3 | -6,-3,0,3,6 |
| 4 | -3/2,0,3/2 |
| 5 | 0 |

translation determinant k(k-3E)(k+3E)の失敗も既に含む。
従って全sourceは

\[
 b\in\{-14,-11,-19/2,-8,-13/2,-5,-2\}
\]

の**七whole a-lines**だけへ還元される。
E零のb=-8をEで割って落とさず、元P0へ直接戻す。
全aのexceptional ranks、a=4又はa=-b-4というnormalizer/leading境界も原証明通り区別する。

## 6. 次の有限判定のscope

全concurrent m2doubleはh=f全空、h=L^2Mの十whole a-lines、h=L^3の七whole a-linesへ尽くされる。
actual linear normalizerの逆、全higher variables、全a wallsを戻した十七族の
ONE有限QQ(a) full30/9/normalizer判定なら、point scanでなく全sourceの有限手順になる。
generic rankだけで全aを覆わず、minor/denominatorの全irreducible wallsを必ず含める。
像24以下、actualDの必要因子、又は新Ore認識の最弱結論で十分。
未閉なら同じprofilesを足さず、保存したwhole-source方程式の構造を変えて読む。

このpureproofと旧全normalform補集合の独立確認前に全G/E完了を宣言しない。
nonreducedR、m3及び他support surfacesは別義務として保持する。
