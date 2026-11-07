# rank零 plane・零symbol全39源の Lie algebroid 次元境界：別作者監査

2026-11-07。reviewer `exceptional_coverage`、producer `plane_collision_review`。
**全段階を独立に読解し PASS とする。** CAS の再実行や数値標本による確認ではなく、
全 source・全補集合を保持した幾何証明の監査である。

受理するのは、rank4 Poisson tensor on P4、実退化五次 D が既約被約、rank零 plane S、
canonical normal-weight-zero symbol が恒等零、という仮定での全射影像の次元上界24。
全 E 帰属、個々の orbit、全 Poisson 分類の完了とはしない。
正確な入力と hash は対応 review job/report に固定した。

## 1. 全39係数からの実 bundle chart を再導出した

既受理の全重複度 conormal 可換性を使うので、旧 m=2 に限定しない。
canonical35/30/9 source の weight零を零とした全30+9を出発点とする。
S=(z=w=0)、base座標をa0,a1,a2とする。
w!=0で t=z/w,xi=ai/w、z!=0で t'=1/t,x'i=xi/t。
よって P4 minus S は Tot(O(1)^3)、dual は A=O(-1)^3。

元の全30係数について商の微分を直接使うと

```
{xi,t}=Pi(ai,z)/w²-z*Pi(ai,w)/w³
      =C_i+tB_i-t²E_i-t³D_i,
{xi,xj}=tLij(x)+Nij(x)
        +xi(t²D_j+tE_j)-xj(t²D_i+tE_i).
```

Π1 の normal-normal 係数は零なので、追加の商微分項は無い。
後者は fiber-linear、Π2 の全九係数は t の二次式を係数とする constant vertical bivector。
anchor は三本の O(3) section、bracket は A の全 differential skew bracket となる。
二 chart の transition は linear Poisson/algebroid の通常の Leibniz transition を持つ。
source を一つの affine chart に切り詰めた還元ではない。

fiber Euler 次数が−1,−2で異なるため、全 Jacobi はそれぞれ
Π1の Lie algebroid 条件と Π2の CE 二cocycle条件になる。
Π2同士の括弧は coefficient が base だけなので自動零。
元 tensor は dense open への restriction で一意に決まり、この写像は単射である。
逆に全 abstract data がP4へ延びることを仮定する必要は無い。

## 2. anchor の全補集合と kernel の型

anchor の共通零 t0 があれば、その fiber A3 上では Π が三次元 vertical 方向だけにあり、
rank≤2。fiber closure は P4 の hyperplane なので actual D に一次因子が生じる。
∞も第二 chart で同じ。恒等零 anchor は global rank≤2であり、いずれも仮定に反する。
従って anchor は全P1で全射、Kの階数二、次数−5となる。

local lift of partial_t の adjoint は K の holomorphic connection であり、Jacobi は
Lie bracket を平行に保つ。小disc上のODEで bracket-constant frame を取れるので、
abelian/nonabelian型は全点で一定である。generic nonabelianだけの議論ではない。

Kが可換なら inner action が零なので、この connection は lift に依存せず大域的となる。
curve上の curvature は零、P1 は simply connected なので K は trivial flat bundleとなり、
deg K=−5に矛盾する。従って全fiberはaff1。

## 3. kernel の全extensionと algebroid の完全復元

aff1の基底を [h,e]=e と取ると、trace character χ(h)=1,χ(e)=0であり、
全 bracket は χ(k1)k2−χ(k2)k1と書ける。
従って L=kerχ=O(-5)、quotient はχで標準的にOと同定される。
全 kernel bundle は H1(O(-5))=C4 の extension class で記述される。
split class は H0K!=0 を与え、K⊂O(-1)^3に反する。
非零 scalar倍は L の真の bundle automorphism で同じ Lie algebra bundleとなり、
全 kernel は P3 の四affine chartで覆われる。未処理の根型は捨てていない。

導分について、D(h)=b*e,D(e)=d*eがaff1の全導分で、いずれもinnerであることを直接確認した。
centerは零。scalar symbolを持つ first-order Lie derivation の bundle Der(K) は
0→adK→Der(K)→TP1→0を満たす。
a→[a,-]|Kは O-linear（Kのanchor零を使う）、JacobiでLie導分となる。
kernel上はad同型、anchor quotient上はidentityなので、A→Der(K)は全同型である。
central extensionなどの別のalgebroid parameterは残らない。

P3の各chartでextensionをČech cocycleとして構成し、基底微分でDer(K)も代数族となる。
全 algebroid は Isom(A,Der(K_e)) によって実際に被覆される。
nonempty fiberは underlying bundle の全 isomorphism の集合で、Aut(A)=GL3のtorsor、次元九。
空fiberを捨てる操作は上界に影響せず、flatness仮定も要らない。
従って全Π1 source dimension≤3+9=12。
baseP1を固定した全eを数えており、normal PGL2の三次元を追加しない。
同時に、frameを商で潰した数でもなく、全GL3を残した source 上界である。

## 4. CE の analytic local-to-global と全九cocycle

local horizontal lift partial_t は bracket-preserving connection を与える。
parallel frameで [horizontal,h]=[horizontal,e]=0となり、rank一baseには曲率項が無い。
従って local algebroid は Tdisc と constant aff1 の直積に同型。
CE complex は holomorphic de Rham と finite CE(aff1) のtensor productである。
aff1では dh*=0,de*=−h* wedge e*なので local CE cohomology sheavesは
q=0,1で constant C、q>=2で零となる。H1の生成元χは全Aut(aff1)で固定される。

cohomology-sheaf spectral sequence の total degree2には
H2(P1,C)=C以外の項が無く、hypercohomology H2は高々一次元。
transgressionがそれを殺す可能性は残すが、強い消滅は要求しない。

ここで differential は O-linear **ではない**。
比較は GAGA が微分写像そのものを比較するとの主張によらない。
各term O,O(1)^3,O(2)^3,O(3) は analytic と algebraic の双方で高次cohomology零であり、
両hypercohomologyは同じ global cochain complexで計算できる。
その global sections はいずれも同じ次数の斉次多項式、微分も同じalgebraic係数の式である。
従って上の analytic H2上界が元の全有限cochain complexに適用される。

global degree1 cochain dimensionは六、degree2は九。
Z2/B2 dimension≤1,dim B2≤6なので、**全九変数の**cocycle空間はdimension≤7。
Π2をshearで全て消せると仮定せず、最大一つの非exact classを保っている。

## 5. 次元復元と受理する範囲

fixed planeで全Π1と全Π2の source は高々12+7=19。
plane choice Gr(3,5) の六次元を戻す。rank4 openのtensorは非零であり、
scalar C*作用はこの全sourceを保つので projectivizationで一を引ける。
従って全像は19+6−1=24次元以下となる。

全anchor共通零、abelian kernel、split extensionの補集合にはそれぞれ明示した矛盾がある。
generic-onlyのparametrization、normal-form標本、39変数GB、higher jetは使っていない。
外部の未読分類定理も追加していない。全重複度conormal補題と旧canonical sourceは
既存の独立受理済み入口への依存であり、今回その全係数行列を再計算したとはしない。
以上の正確な射程で producer の G≤24 命題を受理する。
