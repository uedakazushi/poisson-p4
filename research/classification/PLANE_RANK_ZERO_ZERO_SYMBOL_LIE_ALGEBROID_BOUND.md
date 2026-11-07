# rank 零 plane・weight零symbol零の全源：Lie algebroid による G≤24

2026-11-07。主義務 G、producer `plane_collision_review`。
**命題。** \(\pi\) を \(\mathbb P^4_\mathbb C\) 上の rank 四 Poisson tensor とし、
実退化因子 D が既約被約 quintic、plane \(S=(z=w=0)\) 上でπが零とする。
canonical normal filtration の weight零symbolが恒等零なら、
その全 projective Poisson 像は次元≤24である。
有限標本・normal-coneの一般型・Dのnormalization・G/E帰属を仮定しない。

本稿は純粋な producer 幾何証明で、別作者監査前である。
新しい CAS、39変数GB、jet、normal-form点の抽出は行わない。
full30 weight一源とfull9 weight二源を全て保持し、
全源を rank三 negative bundle 上の transitive Lie algebroid と二cocycleへ写す。
必要なのは dimension上界だけで、Lie algebroidの個々のorbit分類やcocycleの全exact性は要求しない。

## 1. 全39源の実chartとbundle

被約実 D の全重複度 conormal可換性により、既受理の canonical35/30/9 source を使う。
weight零を零とするので全 tensor は \(\Pi=\Pi_1+\Pi_2\)。
base coordinatesを \(a_0,a_1,a_2\)、normal coordinatesをz,wと書くと、全30源は正確に

\[
\begin{split}
 \Pi_1={}&zL+wN+z^2\mathbf D\wedge\partial_w
 +zw(\mathbf B\wedge\partial_z+\mathbf E\wedge\partial_w)
 +w^2\mathbf C\wedge\partial_z.
\end{split}
\]

L,Nはbaseにlinearなbase-base bivectorの全18係数、
\(\mathbf B,\mathbf C,\mathbf D,\mathbf E\) はconstant base vectorの全12係数。
\(\Pi_2\) はbase三pairそれぞれに任意binary quadraticを掛けた全9係数である。
global Euler gaugeの全冗長性を除いた既受理basisであり、未知係数を新たに零としていない。

\(B=\mathbb P^1_{(z:w)}\) とし、\(\mathbb P^4\setminus S\to B\) をnormal projectionとする。
w1 chartで \(t=z/w,x_i=a_i/w\)、z1 chartで \(t'=1/t,x'_i=x_i/t\)。
従ってこのopenは \(\operatorname{Tot}\mathcal E\)、\(\mathcal E=O_B(1)^{\oplus3}\) である。
\(A=\mathcal E^*=O_B(-1)^{\oplus3}\) と置く。
上の全Π1の induced brackets は

\[
 \{x_i,t\}_1=C_i+tB_i-t^2E_i-t^3D_i=:v_i(t), \tag{1}
\]

\[
 \{x_i,x_j\}_1=tL_{ij}(x)+N_{ij}(x)
 +x_i(t^2D_j+tE_j)-x_j(t^2D_i+tE_i). \tag{2}
\]

従ってΠ1はfiberにlinearなPoisson candidate、Π2はfiberにconstantなvertical bivectorである。
両chartのtransitionから、(1)は \(\rho:A\to T_B\) の全三cubic section、
(2)はA上のskew differential bracketを表す。
全Π2は \(C^2=H^0(B,\Lambda^2A^*)=H^0(O(2))^{\oplus3}\)、dimension九。
係数はこの全二chartで変換されるので、normal affine chartだけを使ったsourceの切断ではない。

fiber Euler weightでΠ1,Π2はそれぞれ−一、−二。
Jacobiの三weightは別々だから、全Jacobiは

\[
 [\Pi_1,\Pi_1]=0,\qquad[\Pi_1,\Pi_2]=0,
 \qquad[\Pi_2,\Pi_2]=0
\]

を与える。最後はΠ2がverticalでcoefficientがtだけに依存するため自動零。
第一式はAのLie algebroid構造、第二式はそのtrivial-coefficient Chevalley--Eilenberg differential
\(d_Ac=0\) と正確に同じである。
元tensorはdense openへのrestrictionで一意に決まるので、元の全39源からこのalgebroid/cocycle dataへの写像は単射。
全abstract dataがP4へextendするという余分な仮定はdimension上界に不要である。

## 2. anchorの全零補集合を実Dの因子で除く

anchorが恒等零なら全bivectorは三次元fiber方向だけにあり、rankは高々二。
これは実rank四の仮定に反する。
anchorの三cubicに共通零 \(t_0\in B\) があれば、fiber \(\mathbb A^3\) 全体で
\(\{x_i,t\}=0\)、全Π2もverticalなので、やはりrank≤二である。
このfiberのP4でのclosureは \(t_0\) に対応するhyperplane \(H_{t_0}\) で、fiberはそのdense open。
従って actual quintic Dは \(H_{t_0}\) を成分に持つ。
∞も第二chartにより同じであり、全共通零を扱った。

実Dの既約被約性から、残る全sourceは**全Bでsurjective**なanchorを持つ。
従って

\[
 0\longrightarrow K\longrightarrow A\overset\rho\longrightarrow T_B\longrightarrow0,
 \quad\operatorname{rank}K=2,\quad\deg K=-3-2=-5. \tag{3}
\]

kernel K のbracketはO_B-bilinearなLie algebra bundleである。
local lift of∂tのadjointはK上のholomorphic connectionで、Jacobiによりbracketを保つ。
curve上の小discでparallel frameを取れるので、KのLie algebra型は全点で同じ。
genericだけで非可換とし、退化点を捨てる議論ではない。

## 3. kernelは全fiber aff(1)、その全bundle dataはP3で被覆

Kがabelianなら、\(\nabla_Xk=[a,k]\)、\(\rho(a)=X\) はliftによらず定義される。
Jacobiでflat、holomorphic flat bundleはsimply connectedなP1上trivialなのでdegree零。
(3)のdegree−五と矛盾する。
従って全fiberは唯一の非abelian二次元complex Lie algebra
\(\mathfrak g=\mathfrak{aff}(1)=\langle h,e\mid[h,e]=e\rangle\)。

canonical trace character \(\chi(k)=\operatorname{tr}\operatorname{ad}_k\) は
KからO_Bへの全域surjection。
derived line \(L=\ker\chi=[K,K]\) はdegree−五なので \(L\simeq O_B(-5)\) である。
全bracketは

\[
 [k_1,k_2]=\chi(k_1)k_2-\chi(k_2)k_1 \tag{4}
\]

と、exact sequence \(0\to L\to K\overset\chi\to O\to0\) から一意に決定される。
したがって全Lie algebra bundleはextension class
\(e\in H^1(O(-5))\simeq\mathbb C^4\) で記述される。
e0ならKはOを含み、\(H^0(K)\ne0\) が \(K\subset A=O(-1)^3\) に反する。
Lのscalar automorphismでeの非零scalar倍は同じLie algebra bundleとなる。
全dataは \(\mathbb P H^1(O(-5))=\mathbb P^3\) の有限個の通常affine chartsで被覆される。
binary-cubicの特殊orbitやanchor rankをnormal formの標本へ置き換えていない。

## 4. A→Der(K)の全復元とΠ1 source≤12

\(\operatorname{Der}(K)\) を、K上のfirst-order differential operators Dで
symbolがあるbase vector field Xとなり、
\(D[k,l]=[Dk,l]+[k,Dl]\) を満たすもののbundleと定義する。
local Lie algebra trivializationでは \(\operatorname{Der}(K)\simeq T_B\oplus\operatorname{Der}(\mathfrak g)\)、
anchorはsymbol Xであり、commutatorをbracketとするtransitive Lie algebroid。
\(\mathfrak{aff}(1)\) はcenter零、\(\operatorname{Der}(\mathfrak g)=\operatorname{ad}\mathfrak g\)（dimension二）。
従って

\[
 0\to\operatorname{ad}K\to\operatorname{Der}(K)\to T_B\to0.
\]

\(a\mapsto[k\mapsto[a,k]]\) はJacobiによりLie derivationで、symbolは正確にρ(a)。
Kのanchor零からaにO-linearである。
K上ではcenter零によりad同型、anchor quotient上ではidentity。
(3)との比較で、このmapはbundle及びLie algebroidの同型
\(A\simeq\operatorname{Der}(K)\) になる。
kernel Lie algebra bundleだけから元の全Aを復元するので、別のcentral extension modulusを落としていない。

P3上のextensionは、eの一係数を一にした四affine chartでČech cocycleにより実際に代数族として作れる。
各K_eのLie bracketは(4)、Der(K_e)のtransitionもそのcocycleとbase微分から代数的に作れる。
全Π1はこの族の \(\operatorname{Isom}(A,\operatorname{Der}(K_e))\) で被覆される。
nonempty fiberは \(\operatorname{Aut}(A)=GL_3\) のtorsor、正確にdimension九。
isomorphismの無いeのfiberは空であり、跳躍を無証明のflatnessで無視しない。
finite-type fiber dimension theoremから全Π1 source dimensionは高々

\[
 3+9=12. \tag{5}
\]

ここではbase P1を固定した全eを数えているので、normal PGL2の三次元をさらに加えない。
全frameはIsomの九次元で既に数え、非零tensor scalarもAのhomothetyに含まれる。
sourceを商として小さく数えたものではなく、全源を被覆する実familyの上界である。

## 5. 全二cocycleのuniform上界七：local analytic CEとglobal比較

CE sheaf complexをanalytic topologyのC-vector sheavesとして取る。
anchorは微分を含むのでdifferentialをO-linearとは仮定しない。
小discでanchor∂tのliftを選び、Kのadjoint connectionをholomorphically parallel-trivializeすれば
\(A\simeq T_{\rm disc}\oplus(\mathfrak g\otimes O_{\rm disc})\)。
このCE complexはholomorphic de Rham complex of discとfinite CE complex ofgのtensor product。
holomorphic Poincaré lemmaと直接のg計算により、そのcohomology sheavesは

\[
 \mathcal H^0=\mathbb C_B,\qquad\mathcal H^1=\mathbb C_B,
 \qquad\mathcal H^{q\ge2}=0.
\]

gのCEでは \(dh^*=0\)、\(de^*=-h^*\wedge e^*\) だからH0,H1は一次元、H2零。
H1 generatorはtrace characterでAut(g)が固定する。
またP1はsimply connectedなのでlocal systemとしてもconstantである。

cohomology-sheaf hypercohomology spectral sequence
\(H^p(B,\mathcal H^q)\Rightarrow\mathbb H^{p+q}(B,CE_A)\) を使う。
total degree二に寄与し得るのは \(H^2(\mathbb P^1,\mathbb C)=\mathbb C\) だけ。
\(H^1(\mathbb P^1,\mathbb C)=0\)、\(\mathcal H^2=0\) なので

\[
 \dim\mathbb H^2(B,CE_A)\le1. \tag{6}
\]

ここでd2-transgressionを決定したりH2零と強めたりする必要は無い。
一方、CEの各termは

\[
 O_B,\quad A^*=O(1)^3,\quad\Lambda^2A^*=O(2)^3,\quad\Lambda^3A^*=O(3)
\]

で、全てhigher coherent cohomologyが零。
analytic coherent cohomologyも同じであり、bounded complexのacyclic-term spectral sequenceにより
hypercohomologyは**全global cochain complex**のcohomologyに等しい。
従って(6)は、実際の全九係数cocycle空間Z2とglobal coboundary B2について

\[
 \dim(Z^2/B^2)\le1,\qquad
 \dim B^2\le\dim H^0(A^*)=6,\qquad\dim Z^2\le7. \tag{7}
\]

全Π2 unknownをgeneric fiberだけで数えたものではなく、任意のtransitive Aに共通の上界。
global exact一cochainはfiber translation \(x_i\mapsto x_i+\alpha_i+\beta_i t\)、
すなわちactual linear shear \(a_i\mapsto a_i+\alpha_iw+\beta_i z\) に対応する。
ただし本証明は全cocycleのexact性を仮定せず、最大一つのH2 classを残している。

## 6. 全像 G≤24と明示補集合

(5)と(7)より、fixed planeでの全元Poisson sourceはdimension高々19。
このfamilyは全39係数sourceの実restrictionを含み、anchorの零補集合は§2でactual D factorとして処理した。
plane choiceはGr(3,5)の六次元。
rank四openではtensorは非零で、そのC* scalar方向があるためprojectivizationで一を引く。
全像の上界は

\[
 12+7+6-1=24.
\]

したがって命題が従う。
Lie algebroid moduliの有限点・Π1のgeneric anchor・未知Π2のgeneric rankだけから全源を推論していない。
Π0非零 scalar源、非零 traceless dependent global-V源、他のrankや全Eは本稿で再開・変更しない。
全old packet、global/main/Gitを保持し、新しいpureproofの独立監査へ渡す。
