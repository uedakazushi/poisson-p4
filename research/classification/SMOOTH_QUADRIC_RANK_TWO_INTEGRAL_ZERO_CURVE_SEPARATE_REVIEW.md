# 滑らかな quadric・rank 二・integral 零曲線：別作者の全源監査

2026-11-07。producer `root`、reviewer `plane_collision_review`。
`SMOOTH_QUADRIC_RANK_TWO_INTEGRAL_ZERO_CURVE_EXCLUSION.md` を独立に読解し、
以下の全仮定と全補集合の下で **命題を受理する**。
純粋な幾何監査で、新しい CAS、標本、formal jet 又は行列計算は無い。

**受理する結論。** 実退化因子が既約被約 quintic である rank 四 Poisson tensor について、
非正規曲面成分 S が smooth quadric、\(\pi|_S\) の一般 rank が二、
その anticanonical 零 divisor R が integral reduced となる全枝は空。
R の smoothness、normal-cone の一般 node 性、vertical 不在、Dの正規化又は導手のsmoothnessは仮定しない。
R可約・非被約、quadric cone、他のSは明示残余である。

## 1. prime CI の symbolic power と degree五の全源

\(I=(h,q)\) は smooth quadric の homogeneous prime CI、degree は1,2。
regular sequence の associated graded は
\(\operatorname{gr}_I R\simeq(R/I)[H,Q]\) で、各 graded piece は domain R/I 上 free。
g∉I の乗法は各 piece で単射なので、全 I^m は I-primary。
従って generic S-order≥m は ordinary I^mへの membershipを意味し、symbolic/ordinaryを混同しない。

I3のdegree五では q3（degree六）は使えず、全項がhの倍数となる。
S⊂SingDからm≥2、実D既約からm≥3を排除し、mは正確に二。
I2のdegree五は正確に

\[
 F=h^2A_3+hqB_2+q^2C_1
\]

である。高normal-orderの項はこの全表示に含まれ、後で捨てていない。

## 2. 一般Sへの parallel quadratic の適用

SはPoisson invariant、\(\pi_S\) がunitなdense openでnormal Pfaffian一次係数を見ると
\(\{I,I\}\subset I^2\)。S smoothなのでconormal bracketはbundle morphism、全Sへ零として延びる。
したがって \(I/I^2=N^*\) のHamiltonian connectionはlift非依存、Jacobiでflat。
\(U=S\setminus R\) では \(\sigma=\pi_S\) がsymplecticなのでordinary flat connectionになる。

Hamiltonian flowは \(\pi^2\) と \(\sigma\) を保存する。
normal-order二のassociated gradedとadjunctionを取り

\[
 \gamma=G/\sigma\in\operatorname{Sym}^2N^*\otimes\det N
\]

がparallelとなる。Hamiltonian fieldsはTUを張るので全方向のparallel性。
ここで \(S=\mathbb P^1\times\mathbb P^1\)、\(H=(1,1)\)、
\(\det T_S=O(2H)\)、\(N=O(H)\oplus O(2H)\)、\(\det N=O(3H)\) が整合する。
Rのsingular点を仮定から除く必要はない。
Uはconnected、Gは非零なのでparallel transportによりγは全Uで非零、rank一定。
Gのdivisorial common zeroはRにだけ支えられる。

## 3. primitive root saturation と rank一の parity

関数体上のroot lineをN*にsaturateする。
smooth surface上でsaturated rank一subsheafはreflexive（quotient torsion-freeとdepth lemma）なのでline bundle。
そのdualをAとすればprimitive map \(\ell:N\to A\) を得る。
codimension二のbasepointは許される。

rank一で \(G=\gamma_0\ell^2\) と書くと、各DVRでellの少なくとも一係数がunitだからγ0のpoleは無い。
line bundleのHartogsで全S上regularに延びる。
従って \(V=\operatorname{div}\gamma_0\ge0\)、\(2A+V=5H\)。
parallel性よりVのsupportはR、integral reduced RゆえV=kR、\(R\sim2H\)。
\(2A=(5-2k)H\) はPicS=Z2の両成分が奇数となり不可能。
generic double-root/cuspなどをこのrank一枝から落としていない。

## 4. rank二の分裂と全vertical分岐

\(\gamma\in\operatorname{Sym}^2N^*\otimes\det N\) のdiscriminantはuntwisted scalar。
flatnessからU上非零constant c、従って全sectionとして
\(\Delta(G)=c\sigma^2\)。C上でcは平方なので、関数体の二rootは分裂する。
両rootを上と同じ方法でsaturateしprimitive \(\ell_i:N\to A_i\) を得る。
primitive二linear formの積も各DVRでprimitive（residue field上の非零二linear formの積は非零）なので、
\(G=\gamma_0\ell_1\ell_2\) のscalarγ0は全Sでregular。

\(V=\operatorname{div}\gamma_0\ge0\)、\(A_1+A_2+V=5H\)。
\(\delta=\det(\ell_1,\ell_2)\) は \(A_1+A_2-3H\) のsectionであり、
\((\gamma_0\delta)^2=c\sigma^2\) から
\(\gamma_0\delta=\pm\sqrt c\,\sigma\) という実sectionの恒等式を得る。
従って \(V+\operatorname{div}\delta=R\)。
integral reduced Rの分解なので、全可能性はV0又はV=Rだけ。
vertical成分の無い場合を最初から仮定していない。

## 5. V0：blowup residue と数値Gysinを直接検査

blowup \(b:Y=\operatorname{Bl}_S\mathbb P^4\to\mathbb P^4\)、\(E=\mathbb P(N)\) へπはregularにliftする。
local chart v=uτでのbracketの分子はconormal可換性で割れる。
実退化因子はD'+E、両成分の係数は一であり、inverse formのE-residueβは閉有理一形式。
divisorial poleは \(\Gamma=D'|_E\) だけである。

V0ではΓは二つの相異なるdegree一horizontal prime \(\Gamma_1,\Gamma_2\) だけ。
primitive rootのcodimension二basepointはfiber内のcodimension二集合となり、追加divisorを作らない。
Uでnormal quadraticはrank二なので、generic transverse型はnodeと**結論される**。
generic normal flat frameでnormal inverseは非零係数の \(du/u\wedge dv/v\)。
blowup後のβはfiber座標τについて非零係数のdτ/τを持ち、二primeでのresidueはλ,−λ、λ≠0。
閉性からそれぞれconstantになる。係数体をPoisson中心とは仮定しない。

Gysinの必要な数値部分を、全点SNC又は解消上のpole次数を仮定せず直接確認する。
\(\xi=O_E(1)\)、\(D_0=A_1-A_2\) と書けば
\([\Gamma_1]-[\Gamma_2]=p^*c_1(D_0)\)。
\(A_n=\xi+n p^*H\) をn十分大でvery ampleに取る。
\(L=(1,0),(0,1)\) のそれぞれについて、A_n2と \(A_n(A_n+p^*L)\) の一般CI曲線を選ぶ。
両曲線はΓのsingular locus・二primeのintersection・全bad codimension二集合を避け、Γをtransverseに切る。
compact curve residue theoremより両intersectionに対し
\(\lambda(\Gamma_1-\Gamma_2)\cdot C=0\)。差を取ると

\[
 0=p^*D_0\cdot A_n\cdot p^*L=D_0\cdot L
\]

である（\(p_*\xi=1\)、baseはsurfaceなので三つのpullbackの積は零）。
二つのLでD0の両整数成分が零、従ってA1≡A2。
他方A1+A2=5Hなので2A1=5Hというparity矛盾。
これは既存clean-surface Gysin部分の全scopeと整合し、追加SNC仮定なしに今回必要な数値関係を独立検証する。

## 6. V=R：実全quinticのlinear factor

δは全Sで非零、従って \(N\to A_1\oplus A_2\) はbundle同型。
\(\operatorname{End}(O(H)\oplus O(2H))\) のglobal行列は
\(\left(\begin{smallmatrix}a&0\\c&d\end{smallmatrix}\right)\)、a,d定数、\(c\in H^0(O(H))\)。
rank一idempotentではdiagonalが(1,0)又は(0,1)で、global unipotent conjugationによりcを消せる。
従って全line direct summandはO(H)又はO(2H)で、Pic0 twist等の隠れた別型は無い。

一方のroot mapは \(N\to O(H)\)。
\(\operatorname{Hom}(O(2H),O(H))=0\) より、これは非零constant倍のdhだけであり、Gはh normal coordinateで割れる。
§1の全表示ではq2係数C1のrestrictionが零となる。
Sの唯一のlinear equationはhなのでC1=ch、全actual Fがhの倍数となる。
高normal-orderを捨てた因子推論ではなく、degree五のI2全体に対する因子である。

## 7. 最終判定と補集合

generic rank一、rank二/V0、rank二/V=Rを全て排除した。
parallel tensorの一般smooth Sへの適用、primitive saturation、実sectionのdiscriminant、
全vertical分岐、generic node residue、codimension二を避ける数値Gysin、
global normal summandとactualCI式まで独立に確認し、producer命題を受理する。
R可約・非被約、quadric cone、他次数、全G/E・成分数は閉じない。
global/main/Git、producer証明書とscalar成果は変更していない。
