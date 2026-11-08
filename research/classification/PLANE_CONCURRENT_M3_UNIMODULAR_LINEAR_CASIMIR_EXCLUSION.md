# K0 constant方向・unimodular全源は既約quinticを持たない

2026-11-08、PRIMARY G、producer純証明・独立review前。
新CAS、source probe、matrix、GB又はjetは無い。
全central/全Wを保持し、K0の二fixed族と二period族のt0を排除する。

## 1. 仮定と全対象

P=Q(y)+∂c∧W(y)、Q strictquadratic Poisson、
E_y∧L_WQ=0、actualDは被約既約quinticとする。
Qは先行binary Ore sourceで、U=P2R0,V=Q2R0、R0はconstant normal field。
nonzero linear normalCasimir Lを選びR0(L)=0。
Qの全L-columnは零で、divQ0を仮定する。

先行whole normalformsのB_axis/B_otherはdivQ0。
C_tan/C_transは共に

`divQ=t*(u∂a+v∂b)`

なのでt0の全C2がこの仮定を満たす。
このlemmaをlinearCasimirの無いSS/Jordanへ適用しない。

## 2. ChartL1では全moduleがPoisson vector fieldになる

normal complementをrとしL1で(a,b,r)を用いる。
Q_L-column零なのでZ0、whole fullsourceは

`P_chart=pi+∂c∧Y`, `L_Ypi=0`。

divQ0はこのchartのdiv pi0を与える。
三変数polynomial bivectorのdual一formはclosed、degree≤二なのでpi=Jac(F)となる。
binary sourceのliteral structureから

`F=a Q2(r)−b P2(r)+C3(r)`、degreeF=d≤三。

C3'は元の全central C2で、integration constantはFの定義に吸収出来る。
P2/Q2は同時零でなければよい。
全central sourceをzeroに固定したものではない。

Jacobian bivectorのLie derivative恒等式は

`L_Y Jac(F)=Jac(YF)−(divY)Jac(F)`。

したがって

`d(YF)=(divY)dF`。                                     (1)

## 3. YFはFのpolynomial、generic-fiber仮定を隠さない

例えばQ2非零ならC(a,b,r)=C(F,b,r)という正確なrational coordinate field表示を持つ。
(1)をこのfieldで読むとYFはb,rに依存せず、C(F)へ入る。
Q2零/P2非零ならa,bを交換した同じfield表示を使う。
これでunreadなgeneric-fiber分類を使わずrational functionP(F)を得る。

YFは元chartのpolynomial。
P(t)=A(t)/B(t)をcoprimeに取り、Bezout identityをFへ代入するとA(F),B(F)もcoprime。
B(F)非constantならA(F)/B(F)はpolynomialにならない。
ゆえにBはconstant、YF=P(F)は真のpolynomialである。

## 4. Quintic homogenizationは必ずproper factorを与える

actualDのchart値はsignを除きJac(F) wedge Y、すなわちYF=P(F)。
homogeneous actualDはdegree五なので、F_d=L^d F(a/L,b/L,r/L)をdegree dへhomogenizeして

`D=L^5*P(F_d/L^d)`。                                   (2)

Fがnonconstantならd≥一。
Yはchartでdegree≤三、gradFはdegree≤二だからdegree(YF)≤五。
Pのdegreeは高々floor(5/d)。

* d三ならPはdegree≤一、(2)の全項はL²で割れる。
* d二ならPはdegree≤二、全項はLで割れる。
* d一ならDは二linear forms F_1,Lのbinary quintic。
  複素体上でlinear factorsへ分解し、degree五の既約polynomialにはならない。

P零ならD零でrank四既約quinticの対象外。
Fconstantならpi零、QもL≠0のdense openで零なのでQ零、やはりD零。
各zero complementを含め、全sourceは既約quinticを持たない。

## 5. 結論と残余

K0 constant-directionのB_axis/B_other全fixed source、
C_tan/C_transのt0全magnetic sourceは**空**。
全actualW/central coefficientsを保持した結論である。
C_tan/C_trans t非零、SS/Jordanのprimitive invertible linear normal方向はまだ残る。
SSのunimodularlineやJordanのunimodularpointへ、linearCasimirを未確認のまま本lemmaを拡張しない。
全K0、concurrent G全体、E認識又はcomponent数はまだ未完了。
