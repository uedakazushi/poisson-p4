# cusp接線・割線の全lambda：相異なる22次元軌道と厳密なGa二乗安定化群

2026-11-08、PRIMARY E、`/root/log_manuscript`。純係数証明。
既固定source／actual normalizerを使い、新operator、point、rank、minor、profile、jet又はCASは無し。
新component又は既知親への陽性包含を主張しない。

## 1. 全normal formと固定Dの射程

全六磁気sourceのactual有限shearは
`TWO_LINES_CUSP_TANGENT_SECANT_MAGNETIC_NORMAL_FORM.md` に固定される。
磁気I単元では、座標(h,x,y,z,w)で

\[
 P_\lambda=U\wedge\partial_z+V\wedge\partial_w+
       3(xz+\lambda h^2+y^2)\partial_z\wedge\partial_w,
\]

\[
 U=-hx\partial_x+(3(h-x)^2-3hy)\partial_y,\quad
 V=3x(h-x)\partial_x-9xy\partial_y
\]

の全lambda線へ戻る。
actual Dは全lambdaで \(-9hxC_3\)、\(C_3=h^2y-(h-x)^3\)。
以下は、この全lambda源の**同値関係**を直接解く。
任意の同値行列のblock形はDのessential kernelから導くので、
固定Dを持つ未知sourceの分類又は全ambient親を追加仮定する必要は無い。

## 2. marked base配置の射影自己同型は自明

\((X,Y,Z)=(h-x,h,y)\) とすれば、cubicは標準cuspidal cubic
\(Y^2Z-X^3\)、二直線は \(Y=0\) と \(Y-X=0\) となる。
前者はcusp接線、後者はcuspを通る非flex割線であり、contact型が違うので交換されない。
唯一のdegree三素因子も保たれる。

cubicのnormalizationを \([t^2:t^3:1]\) とする。
cuspはt零、唯一のsmooth flexはt無限遠であり、全plane automorphismはその両点を固定する。
従ってnormalization上ではt→s t、s単元であり、そのambient変換は
\((X,Y,Z)\mapsto(s^2X,s^3Y,Z)\)。
この変換がY-Xを同じlineへ送るにはs³=s²、すなわちs一でなければならない。
cubicをpointwise固定するplane projectivityは、同曲線がP2を張るのでidentityである。
よって \(\operatorname{Aut}(\mathbb P^2,hxC_3)=\{1\}\)。
これは有限という上界だけでなく、全base marked配置の自明性である。

## 3. 全PGL同値は全block変換と磁気cokernelで決まる

PλとPμがactual PGL5で同値なら、その退化五次も同値である。
Dのdirectional derivative kernelはfiber(z,w)の二次元空間そのものでintrinsic。
従って全PGL同値はそのkernelを保ち、base quotient上で§2の自己同型を与える。
global scalar matrixを除けば、実際の全行列は必ず

\[
 b'=b,\qquad f'=Cf+Lb,\quad b=(h,x,y),\quad f=(z,w),\quad C\in GL_2
\]

の形である。baseとfiberの非表示の混合方向を落としていない。
両source及びこのpushforwardはh-Casimir代表なので、同じprojective bivectorである条件は
全raw tensorの非零scalar倍一致と同値である。
これはh一chartがprojective tensorを忠実に表すことからも従う。

\(g_*P_\lambda=sP_\mu\) とする。
mixed行のSaito二場U,Vは定数上で独立なので、全mixed係数の比較から
\(C=sI_2\) が強制される。
したがって残る全actual mapは

\[
 z'=sz+\ell(h,x,y),\qquad w'=sw+k(h,x,y)
\]

だけである。ell,kは任意一次式。
base-base括弧が零なので、有限pushをそのまま計算すると、tensorをsで割った後の全fiber bracketは

\[
 3xz'+3sR+\delta(\ell,k),\quad
 \delta(\ell,k)=-3x\ell+V(\ell)-U(k).
\]

これはinfinitesimal像だけではなく全有限式である。
ell=a h+b x+c y、k=d h+e x+f yを代入した全像は

\[
 \delta=-3f h^2+(-3a+3b+e+6f)hx+3fhy
             +(-6b-3f)x^2-12cxy. \tag{1}
\]

R=r0h²+r1hx+r2hy+r3x²+r4xy+r5y²の全二cokernelは
\(I=r_5\)、\(J=r_0+r_2\) である。
式(1)によるshiftは両方を不変にし、全scalar-fiber mapは
\((I,J)\mapsto(sI,sJ)\) と作用する。
Pλ／Pμでは(I,J)=(1,lambda)/(1,mu)なのでs一、mu=lambdaを強制する。

従って**相異なる有限lambdaのPλはactual PGL5同値ではない**。
canonical sourceの磁気parameterの違いをtensor scale又は非線形座標変換で同一視していない。

## 4. 全lambdaの安定化群と軌道次元

Pλを固定するには§3からs一、式(1)零が必要十分。
順にf零、b零、c零、e=3aとなり、a,dだけが自由である。
よって全安定化群は厳密に

\[
 (h,x,y,z,w)\longmapsto(h,x,y,z+a h,w+3a x+d h),\quad(a,d)\in\mathbb A^2.
\]

二generatorは \(h\partial_z+3x\partial_w\) と \(h\partial_w\)、可換である。
逆は(a,d)の符号反転、det一、全lambdaに同じ式である。
base marked配置に非自明な有限自己同型が無いので、未表示の有限extensionも無い。
従って全有限lambdaでstabilizerはGa²、actual orbit次元は正確に24−2=22。
lambda=-1のzero planeの出現も、この全stabilizer比較を変えない。

全fixed-D sourceとD軌道の次元を合わせた既証明により、全A6の飽和閉像は23次元。
この稠密開のreduced orbit商は、一つのparameter lambdaを持つ。
全I零boundaryの商、componenthood又はscheme smoothnessをこの主張から推測しない。

## 5. 次のactual parent曲線のための停止条件

非自明なone-parameter subgroupの、**固定した一tensorからのprojective極限**は、
そのsubgroupに固定される。従ってそのendpointのPGL stabilizerはGmを含む。
一方、上の全Pλのstabilizer Ga²には非自明なalgebraic torusが無い。
よって固定tensor及び固定conjugationのone-parameter subgroup極限だけでは、
どの有限lambdaのPλも作れない。

これは全known親又はmoving sourceからのactual曲線を否定しない。
period／smoothing parameterやregularなpreconjugationが曲線に沿って動く場合は別である。
したがってこのlemmaを固定frameやpointを変えて試す新profile列の根拠にはしない。
全lambdaのT37陽性曲線、又は全prior比較を伴う新成分性が残る正確なE義務である。
