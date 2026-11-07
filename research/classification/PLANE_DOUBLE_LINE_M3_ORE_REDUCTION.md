# m三squarefree残二源：全Ore sourceを四canonical moduleへ

2026-11-07、主義務G。producer手証明。
Cloud resetで失われた初回19-source receiptはhash138f21ee…として別に記録し、fresh REPLAY receiptで全sourceを再確認する。
本noteは次のfinite module計算を指定するもので、そのrank成功はまだ主張しない。

## 1. 全kernelからliteral Euler gaugeでc-dependenceを除く

残るfirst symbolsは、base(a,b,c)、normal(z,w)で

`P0_D=∂c∧(a²∂a−2ab∂b−2az∂z−2aw∂w)`、

`P0_N=P0_D+bz∂c∧∂w`。

前者の完全weight一kernelは二十、後者は十五。
両方ともP1のab bracketはc-independentであり、c-dependenceはac/bc bracketの
`c A(z,w)`、`c B(z,w)`だけである。
全P2九係数もc-independent。
L=−A∂a−B∂bとしてEuler_total∧Lを引くと、同じprojective tensorを正確に

`Π=Q(a,b,z,w)+∂c∧W(a,b,z,w)`

で表せる。全source係数間のlinear inverseも以下に含める。
Qに全P2のab三係数f2を保持し、残るac/bc六係数もWへ全て戻す。

全projective Jacobiをc次数で比較すると、必要十分条件は

`[Q,Q]=0, E_y∧[Q,W]=0`, y=(a,b,z,w)

である。c-independent strict liftをglobal divergenceから仮定するshortcutではない。
二十九又は二十四の全higher sourceをliteral Euler gaugeで移した結果である。

## 2. D源の全Qとstrict Jacobi

A=A0z+A1wを任意linear form、X=(Xz,Xw)を全quadratic normal vector fieldとして、全Qは

`Qab=−bA+f2`, `Qan=−A E_n`, `Qbn=X`, `Qzw=0`

となる。A二係数、X六係数、f2三係数が全十一Q変数。
完全kernelの残十二係数とfullP2のac/bc六係数がWの全十八変数である。
joint coefficient mapはconstant invertibleであり、sliceを選んで全supportを推論しない。

strict Q Jacobiは正確に`X(A)=0`というbinary-quadratic三係数条件。
実際`[−A E_n,X]=−A X+X(A)E_n`、Qabのb微分は−Aであり、全f2はこのJacobiへ寄与しない。
他の二normalを含むtriplesも零なので、条件は必要十分である。
leading actual termは符号を除き`G3=3a² det(n,X)`。

A零ならQ=∂b∧(X−f2∂a)、L_∂bQ零。
全Πはc-independentで、既constant-factor suspension lemmaによりactual Dはcとbの二方向に依存せずessential-three。
このwhole枝は既八23成分に入る。higher fiberを小さいと仮定しない。

A非零ならactual constant normal GL2でA=zへ直す。
X(A)=0からX=R2(z,w)∂w。
squarefree G3はR2がsquarefreeかつzを因子に持たないことと同値。
normal Borel fixing zは他二rootをw=0及びw=zへ送り、actual b-scaleで共通振幅を一へ戻せる。
P0_Dはnormal GL2及びこのb-scaleで保たれる。
従って全open sourceで

`Qab=−bz+f2, Qaz=−z², Qaw=−zw, Qbz=0, Qbw=w(w−z), Qzw=0`。

## 3. D源の全central extensionは二fixedと一wholeline

actual a'=a+u(n)、b'=b+v(n)でu,vを全linearとすると、

`f2→f2−X(u)`。

vの項は−zv+zE_n(v)で相殺する。
従ってcoboundaryの全imageはC*R2で、残る全classはαz²+βzw。
両係数を一律0/1へ固定しない。

actual diagonalとtensor scalar

`(a',b',c',z',w')=(s a,s b,r c,(s²/r)z,(s²/r)w)`, `σ=s/r`

はP0_D及びcanonical anchorsを保ち、f2を(r/s)f2へ送る。
全residual sourceは、zero、α零/β非零、α非零の三whole strata、すなわち

`f2=0`, `f2=zw`, `f2=z²+τzw (τ∈C)`

で尽くされる。finite root permutationsによるτの重複は許し、whole complex lineを残す。
base/normal shearsで生じるWのEuler二冗長性は`W→W+ell(n)E_y`で全十八gaugeへ戻す。
これはQを変えない正当な五変数Euler gaugeである。

Qの全parameter dimensionは、zero central classで六、nonzero fixed classで七、wholeτで八以下。
A二、R2三、central coboundary一、nonzero central amplitude一、residual modulus一から得た安全なupperである。
W fibersはactual normalizerで同型なので、fixed canonical moduleのfiber dimensionをkとすると、
P0_D source五、plane六、projectivization一を戻した全像upperは

`16+k`（zero）、`17+k`（zw）、`18+k`（wholeτ）、`17+k`（finiteτ wall）となる。

## 4. N源の全Qとsquarefree open

全Qの非共通係数をu,v,r、全central quadraticをf2とする。

`Qab=buz+f2`、

`Xa=(2uz²,3uzw−vz²/6)`、

`Xb=(vz²−6uzw,rz²+2vzw−9uw²)`、`Qzw=0`。

Qの全六係数（u,v,rとf2三）及びW十八が元の全二十四higher係数へ可逆に対応する。
直接計算すると`[Xa,Xb]=uz Xb`なのでstrict Q Jacobiは**全u,v,r,f2で恒等零**。
PfQは`−(v²+12ur)z⁴/6`である。

T=3aw+bzとすると、符号を除いた全leading cubicは

`G3=z[3ra²z²+va zT−uT²]`。

そのdiscriminantはnonzero constantを除き`u²(v²+12ur)a⁸`。
従って指定squarefree openは正確に`u(v²+12ur)≠0`。
零補集合はtarget leadingを持てないので、無根拠にgenericだけを残したlocalizationではない。

## 5. N源は一canonical Q、central f2も全て除去

actual normal shear w'=w+h zはP0_Nを保ち、

`u'=u, v'=v+6hu, r'=r−hv−3h²u`。

Δ=v²+12urは不変。
uΔ非零でh=−v/(6u)としv零、r'=Δ/(12u)非零へ直す。
actual diagonal及びtensor scalar

`(a',b',c',z',w')=(s a,(s l/k)b,r0 c,kz,lw)`, `σ=s/r0`

はP0_Nを保ち、u'=s²u/(r0k)、r'=s²l²r/(r0k³)。
nonzero scalesと有限square-root choiceでu=r=1へ直せる。
全inverseはminus shearとnonzero diagonal inverseであり、higher sourceへ同時に作用する。

canonical anchorsは

`Xa=(2z²,3zw), Xb=(−6zw,z²−9w²), Qab=bz+f2`。

a'=a+αz+βw、b'=b+γz+δwによる全central変化は

`(γ−β)z²+(2δ+6α)zw+9βw²`。

そのimageはSym²C²全体であり、全f2をzeroへ移せる。
したがってN源のsquarefree全支持は一fixed Qで覆われる。
全Q-source dimensionはu,v,r三＋f2三の六。
P0_N source七、plane/scalar五を戻し、全像upperは`18+k`。

## 6. 全W moduleと最弱のfinite target

全WはW0（P0_D又はP0_Nの括弧内field）に次を加えた十八変数である。
Wa,Wbの各々はbase a,bにlinear・normalにlinearな四係数と、全normal quadratic三係数。
Wzはzw,w²二係数、Wwはz²,zw二係数。
raw positive W二十からEuler二方向を固定した完全源で、全P2を保持する。

四canonical Qについて全`E_y∧[Q,W]=0`は140×18 affine linear system。
全original350との同値、complete kernel/minor、particular又はcokernel contradiction、
全actualD affine coefficientsを保持する。
Q lineではmatrixとaugmented matrixのgeneric minor両方の全irreducible factorsを尽くし、
exact characteristic-zero fieldsで全finite wallsを処理する。

十分なfiber targetsはD zeroでk≤8、D zwでk≤7、D genericτでk≤6、finiteτwallでk≤7、N fixedでk≤6。
あるいはaffine actualDの全coefficientsにproper factorがある、全Dがessential-three、又はwhole module不整合ならそれで十分。
未達なら全actual affine familyを保存し、generic rankやfinite-sliceだけで閉じない。
新しいjet、minor列又は大きいGBはこの有限全source targetには不要。
