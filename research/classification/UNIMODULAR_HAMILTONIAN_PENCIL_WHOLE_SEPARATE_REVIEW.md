# unimodular全linear-Hamiltonian／regularpencil／零factor入口：独立査読

2026-11-08、integrator `/root`。旧Gを再開せず、新Eの全低次fiber被覆を手で確認した。新operator/点/minor/jet無し。全既45係数kernelと全Segre24+6係数は旧独立exact証明書を入力にし、本runの再計算とは呼ばない。

## 1. 全linear-Hamiltonian源

元P=C+L+Q、divP0、PfP=kunit、Q#dell0。全affineHamiltonian(V,ell)の零点無し/trace0/classificationとnil六対、全45係数の正確な源を旧noteで照合した。非nilsemisimpleはV=a(x∂x-y∂y)+∂w、semisimpleとconstant作用を有限二次tensor表現で分離し、全weight源にdiv/Pfを課すとproductfamilyへ入る。
constantVは全A/BcrossがquadraticHamiltonianH≤3からなり、literal P_e=(a+ex²)Jxy+bJzw+∂y∧X_H は全HでPoisson。係数のx依存を保ち、全eunitがrawCA全源、閉性でe0全CA。
別々の二planeのquadraticf,gについてP=A J1+B J2+K Xf∧Xgは全parameterでPoisson/PfAB。両Hessianunitのdense開でactualaffine normalizerを用い、全五係数(xy=p xy+A,zw=q zw+B,xz=kxz,xw=-kxw,yz=-kyz,yw=kyw)へ。真T16 period/edge正逆は(q+2k,-p-2k,p-2k,A,B)と逆p=(rho2-rho1)/2,q=rho0+(rho1+rho2)/2,k=-(rho1+rho2)/4。全開pqunitでnamedsource、pq0も同じpolynomialfamilyの閉性。任意quadraticf/gを全体既約源のdense開から戻せる。HQQ全分類は不要。

nil五非空型についてvolumeconstantな補助(t,r,s,l)でV=∂t、全P'=−∂t∧∂l+k0∂r∧∂s+∂t∧Y、divrsY0。この全YにB'=−l∂t∧∂l+∂t∧D、divrsD0がuniversalcompatible、divB'=∂t、Pf(P'+eB')=PfP'*(1+el)。補助でYが高次になっても切り捨てない。元座標へ戻したliteralBがquadraticであることがactual曲線所属の根拠。
特にnil22ではt=y,r=x-y²/2,s=z-wy,l=w、V=(y,1,w,0)、∂l=(0,0,y,1)、D=2(r∂r-s∂s)/3。−wV∧∂lのxz=-wy²とV∧Dのxz=-2(yz+xw)/3+wy²が相殺する。yzは−wy−2z/3+2wy/3=-(2z+yw)/3、xy=(y²-2x)/3。他三括弧は(-yw,-w,-w²)。全六がproducer表に一致。
nil21short ell=yではV=y∂x+∂z、∂l=z∂x+∂y、−yV∧∂l=(-y²,yz,0,y,0,0)。他三行はD0で直ちに同表。全補助Jacobian±1で、div/Pf変換を正当に使える。
全eunitのmodularは正確にeV、actualD=h^4(h+eell)でmarkedh multiple、元PfQ0仮定を満たす。nil31longの一次Pfは旧全A源T01、nil22の一次Pfは旧全A−1/22成分。後三行はrank一translation全known有限閉和（より具体的M1）へ。nil31shortは旧全45係数の空性。全e0元Pも同じ既知閉和へ。simpleHの同rank fieldへこの定理を拡大しない。

## 2. 全primitive regular二次pencil

旧全Segre十型を照合した。全beta=df∧v+u∧dgの8parameterとOmega=J-u∧vによりalpha=Omega+dF∧dG。旧rank2Omega補題はf非退化だけを使い、Omega∧u∧v0からJdegenerate矛盾。八型の全OmegaPf恒等零と四相異固有値空も保持。残SS二block/nil[2,2]のみ。
Omegaunit,Hessfunitでx0=−Hessf^-1*uのactualtranslationを使うとFlinearが消える。{F,G}_Omega constantのlinear項が(Omega^-1 Hessf)x(v'.x)0を与え、可逆性からv'0。全u/vを除算するchartを選ばず同時除去出来る。
SS二blockの逆tensorは上のproductT16。nil[2,2]の全Omega=t(dx1∧dx4+dx2∧dx3)+s dx2∧dx4,tunitについてA=(K,-sK/t;0,K),B=(0,K;0,0),K²=-I/t²。block演算A^n=(K^n,nK^(n-1)(-sK/t);0,K^n)からA³+A/t²=(2s/t³)Bを独立確認した。
sunitでodd polynomial p(z)=t³(z³+z/t²)/(2s)。全T∈sp(Omega)についてOmegaTとOmegap(T)は対称、Tとp(T)可換。F_T=x^T OmegaT x/2,G_T=x^T Omegap(T)x/2、alpha_T=Omega+dF_T∧dG_T は閉でPfconstantunit。inverseは4次元skew complementarymatrix/PfだからliteralquadraticP_Tでparameterに多項式。sp(Omega)のdense regularsemisimple開で±a,±b,ab(a²-b²)unit、quadraticF/Gは二productのlinearcombinationとなりT16。全Tの閉性でAへ戻す。s0はOmega_s全polynomialinverse族の閉性、元u/vはactualtranslation逆で全復帰。
非退化member無し、gcd非constant、gamma0の全補集合をgenericpencilの平滑化で覆ったとはしない。

## 3. 零最高とfactor最高の全入口

gamma0ならJ∧beta=beta²0。betaのimageはJprimitive五次元quadraticspaceのtotalisotropicでdim≤2。dim≤2のdecomposable二形式線形空間は共通constant1form aを持つ。beta=a∧b(x)、i_vJ=a、a(v)=0、J∧a∧b0からb(v)0、ゆえにi_v alpha=a constant。inversePはconstantlinearHamiltonianを持ち全CAへ。beta0も保持。
共通factorのあるgammaがirreducquadfactorだけならgamma=q Omega0、gamma²0でOmega0=a∧b、closednessでqは二変数binaryquadとなりC上linear分裂矛盾。従ってgamma=ell eta。etaimageは全decomposable/wedgezero線形空間、Grassmannian二rulingから共通a因子又はconstant三次元W支持。第一はvolume-dualQのlinearCasimirで上の全membership。第二はi_vgamma0,dgamma0からLie_vgamma0で全三変数pullback。
後者fullpullback及びnonemptylogpoleの全低次延長はまだ未解決。Gのprimitive像排除をE帰属へ代用せず、六最高閉包の全逆像へのroutingのみ受理。

新成分増分零。以上のlinear-Hamiltonian全source、primitive regularpencil全source、gamma0全sourceとfactor全routingを受理する。全unimodularh5/全E/総成分上界は未完。
