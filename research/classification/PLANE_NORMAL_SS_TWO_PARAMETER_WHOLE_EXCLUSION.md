# rank零plane・normal半単純の三つの全二parameter源は実既約quinticでは空

2026-11-08、PRIMARY E。純幾何 producer proof。新 M1／D2／kernel／minor／GB／profile は実行しない。
既存全 kernel と元350行の compatibility から、全係数を保った actual factor 排除を行う。

**定理。** rank零planeに沿う actual 被約既約quinticの位数が二の場合、
normal半単純の全 source \(\mathrm{SS\_MIXH},\mathrm{SS\_MIXN},\mathrm{SS\_NILREG}\) は空である。
全 parameter／finite walls／shear stabilizer jump と inhomogeneous \(P_2\) を保持する。
三parameter SS_MIX、nilpotent normal、low21、他の co-Higgs source は本稿の明示補集合。

## 1. 全 first symbol と保存済み有限被覆

旧 `PLANE_RANK_ZERO_GLOBAL_VECTOR_SEMISIMPLE_FIRST_JETS.md` と
`PLANE_RANK_ZERO_M2_COMPLETION_SEPARATE_REVIEW.md` の全source theorem を使う。
\(H=\operatorname{diag}(1,1,-2)\)、\(N=b\partial_a\)、\(J=b\partial_a+c\partial_b\) とし、
全二parameterの実表示は
\[
 P_0=T\wedge E_n+V\wedge(z\partial_z-w\partial_w)+kzw\partial_z\wedge\partial_w,
\]
\[
 \begin{array}{c|cc}
 &T&V\\\hline
 \mathrm{MIXH}(t,k)&tH+N&H\\
 \mathrm{MIXN}(t,k)&H+tN&N\\
 \mathrm{NILREG}(t,k)&tJ+J^2&J
 \end{array}
\]
全 \((t,k)\in\mathbb C^2\) を保つ。
三者とも \(T,V\) は独立で、\(P_0\) の actual Pf は非零定数倍の
\(zwb^2c\)（mixed）又は \(zwc^3\)（nilreg）である。
従って higher が actual shear だけなら被約既約 D を与えない。

全 \(P_1\) 三十、\(P_2\) 九、全350行の \(D_1,D_2\)、actual shear六列を保持する。
旧 TWOPARAMETER／NIL_LINE reports と別作者reviewは、全 generic・whole line・全根を exact に覆う。
本稿ではその complete kernels／rank certificates を再計算せず使用する。

## 2. normalSS全源の actual shear stabilizer

w一 chart、\(t_0=z/w,x=(a,b,c)/w\) と書き、source parameter t と区別する。
\[
 \Pi_0=E_x\wedge U_x+t_0(Bx)\wedge\partial_{t_0},
 \quad U=T-V,\quad B=2V+kI.
\]
actual translation \(x\mapsto x+v+t_0u\) の一次像が零となる条件は
\[
 Bv=Bu=0,\qquad
 v\wedge Ux+x\wedge Uv=0,\qquad
 u\wedge(U+B)x+x\wedge Uu=0.                                  \tag{1}
\]
有限作用は homogeneous linear shear \((a,b,c)\mapsto(a,b,c)+wv+zu\)、逆は符号反転。

trace零matrix \(R\ne0\) について \(v\wedge Rx+x\wedge Rv=0\) は、
\(R\) が rank 一 nilpotent で \(v\in\operatorname{im}R\) の場合だけ。
これは \(x=v\) から \(Rv=\mu v\)、さらに \(R-\mu I\) の像が \(\mathbb Cv\)、trace零から \(\mu=0\) として従う。
(1) の u式は \(Bu=0\) の下で \(T+V\) の同じ式へ帰着する。
従って二parameter mixedでは shear rankが五になるのは \(t=\pm1,k=-2\) だけ、
nilregでは \(t=\pm1,k=0\) だけ。
MIXNは全域で rank六。全rank低下を明示し、無料の六shearを例外点で引かない。

一般に actual shear rankが \(r\)、\(\dim\ker D_1=r\)、\(\ker D_2=0\) なら、
有限shearで全 \(P_1\) を零にし、変換後の全 \(P_2\) は \(D_2P_2=0\) から零。
有限shearが作る二次補正を捨てていない。
従ってその全 fiber は \(P_0\) と同型で、上記 reducible／nonreduced Pf のため本対象では空となる。

## 3. MIXNの全二parameter

旧全九lineと全rootsで、この剛性を使えない点は正確に
\((t,k)=(3,3),(-3,-3)\)。いずれも complete kernel dimension七、shear六、\(D_2\) injective。
正の点の旧 \(D_1\) kernel最終columnは、旧三十basis row二十九のみ一、
\[
 R=w^2\partial_c\wedge\partial_z.
\]
その mixed anchorは t零で非零、shearの \(t_0B\theta\) は \(t_0\) で割れるため、Rはshear像の外にある。
従って全 \(P_1\) は有限shear後 \(\rho R\)、\(\rho\in\mathbb C\)。
Rは strict Poissonなので weight二forcingは零、injective \(D_2\) より \(P_2=0\)。
この全 sourceは hyperplane b零を保つ。\(P_0\) のb-legはbで割れ、Rはbを動かさないためである。
actual Pf はb因子又は零、矛盾。
normal swap \(z\leftrightarrow w\) と base \(N\mapsto-N\) は \((t,k)\mapsto(-t,-k)\) を実現するので負の点も全逆付きで覆う。
従って MIXN は全域で空。

## 4. MIXHの全十二非剛性点

旧全十line／全rootsの complete certificates を実 parameterに戻し重複を除くと、
非剛性点は次の十二点だけである。その他の全fiberは (24,9,6) の rank certificateと第2節で空。

|点（t,k）|shear後の全残余|factor|
|---|---|---|
|\((1,4)\)|\(P_1=\rho w(2a\partial_a+b\partial_b)\wedge\partial_c, P_2=0\)|b|
|\((-1,4)\)|上のwをzへ替えたもの|b|
|\((3,-2)\)|\(P_1=\rho w^2\partial_a\wedge\partial_z,P_2=0\)|b|
|\((-3,-2)\)|\(P_1=\rho z^2\partial_a\wedge\partial_w,P_2=0\)|b|
|\((1,-8)\)|\(P_1=\rho cw\partial_a\wedge\partial_b,P_2=0\)|c|
|\((-1,-8)\)|wをzへ替える|c|
|\((1,1),(-1,1),(0,1)\)|\(P_1=0,P_2=\rho\{w^2,z^2,zw\}\partial_a\wedge\partial_c\)|b|
|\((0,-2)\)|\(P_1=0,P_2=\rho zw\partial_a\wedge\partial_b\)|c|
|\((1,-2),(-1,-2)\)|全三変数 \(P_1\) と全九 \(P_2\)：第5節|c|

表の前六点は ker dimension七／D2 injective、shear六。
表示Rは全て strict Poissonなので、weight二 equationは \(D_2P_2=0\)、従って \(P_2=0\)。
complete kernelのextra方向であることも復元する。
\((1,4)\) では \(B=\operatorname{diag}(6,6,0)\)。cの二shearは旧column零／一を張り、
他の四shearは四mixed anchorを張る。表示Rはcolumn二でこれらから独立。
\((3,-2)\) のRはBの像に無いa-anchorを、\((1,-8)\) のRは可逆Bの下で零anchorの非零方向を持つ。
負tの点は normal swap、tensor scale−一、base \(N\mapsto-N\) の actual逆付き変換で戻る。
表の四つの D2核一点sourceはP1全六shear、保存された一columnの全P2をそのまま保持した。
各factor欄の hyperplane は full tensor に invariant。
rho零を落としたり、一つの quadratic coefficientだけを選んだ検査ではない。

## 5. MIXH critical点：全 inhomogeneous P二の直接排除

\((t,k)=(1,-2)\) では \(D_1\) kernel八、actual shear五、D2 kernel一である。
旧kernel columnsを零から \(k_0,\ldots,k_7\) とすれば、shear像は
\[
 \operatorname{span}(k_0,k_1,2k_2+k_3,k_6,k_7).
\]
これは (1) からも五 independent actual shearとして確認出来る。
全残余三方向は \(k_2,k_4,k_5\) と取れ、full sourceは finite shear後
\[
 P_1=w(rb-2sa)\partial_a\wedge\partial_b
        +scw\partial_b\wedge\partial_c
        +u zw\partial_a\wedge\partial_z,
 \qquad r,s,u\in\mathbb C.
\]
\(P_2\) を kernel一だけに最初から制限しない。全九係数を
\[
 P_2=A_2(z,w)\partial_a\wedge\partial_b+
      B_2(z,w)\partial_a\wedge\partial_c+C_2(z,w)\partial_b\wedge\partial_c
\]
と保つ。
このとき元の **projective** Jacobi四vectorの二成分は、c零で literal に
\[
 (E\wedge\operatorname{Jac}(P))_{a,c,z,w}|_{c=0}=6zw B_2,
 \qquad
 (E\wedge\operatorname{Jac}(P))_{b,c,z,w}|_{c=0}=6zw C_2.          \tag{2}
\]
手展開を示す。raw cyclic JacobiをJと書くと
\[
 J_{acz}|_{c=0}=z(-2B_2+C_2-2w\partial_wB_2),\quad
 J_{acw}|_{c=0}=w(C_2+2z\partial_zB_2),\quad J_{czw}=0,
\]
\[
 J_{bcz}|_{c=0}=-2zC_2-2zw\partial_wC_2,\qquad
 J_{bcw}|_{c=0}=2zw\partial_zC_2.
\]
Euler degree二から (2) を得る。
strict Jacobi成分を勝手に零としたのではなく、元350行のprojective equationsの必要な組を使った。
polynomial域と非零定数六から \(B_2=C_2=0\)。全r,s,uとA2を保持して c零はinvariantとなり、
actual Pfはcで割れる。これが critical全sourceの排除である。
負のcritical点は前節の実変換とその逆で覆う。

## 6. NILREGの全K零線と全例外

\(k\ne0\) は旧全line／root certificatesで全 (24,9,6)、第2節で空。
\(k=0\) whole t線には旧 NIL_LINE report の全 actual slice／元cokernel350式を使う。
genericではshear六、slice二変数 \(u_0,u_1\)。
旧二つのsquare compatibilityから、全 reduced supportで \(u_0=u_1=0\)。
basis/minor/resultantの全例外は正確に \(t=\pm3,\pm1,\pm1/3\) で、以下で全て戻す。

全tの D2 kernelは \(Q_2(z,w)\partial_a\wedge\partial_b\) の全三方向。
genericと \(t=\pm1/3\) では二つのslice squaresが全P1残余を零にするため、
\(P_1=0\)、全P2がこのkernelとなり c零は invariant。

\(t=3\) では旧full compatibility rows132／134は
\[
 u_0^2/4,\qquad u_0u_2/2+3u_1^2/4.
\]
従って全complex supportで \(u_0=u_1=0\)、\(u_2\) は全て保持される。
残る旧slice columnは row二十七のみ一、
\(P_1=u_2 w^2\partial_a\wedge\partial_z\)。strict self-Jacobiは零、
全inhomogeneous weight二 equationは従って \(D_2P_2=0\) へ戻り、全三方向Q2を保持する。
いずれのtensorも c零を保つ。

\(t=1\) ではshear rankが**五**、slice四変数である。
旧full compatibility row33の \(u_3^2/2\) と row132の \(u_0^2\) は \(u_3=u_0=0\) を強制する。
残りの全二変数の first weight は旧sliceそのまま
\[
 P_1=w[(u_1c+u_2b)\partial_a\wedge\partial_b+u_2c\partial_a\wedge\partial_c].
\]
これは aに係数が依らない \(\partial_a\)-wedgeでstrict Poisson。
従って全P2は再びQ2 kernelで、c零が invariantとなる。
追加のshear stabilizer一方向を無料に引いたり、P2のparticular forcingを無断で零にしたのではない。
生き残った全P1の self-Jacobiが零だからforcingが消える。

normal swap \(z\leftrightarrow w\) と base \(J\mapsto-J\)、\(J^2\mapsto J^2\) は
\((t,k)\mapsto(-t,-k)\) の実可逆変換なので、負の三例外も全て覆う。
以上により全K零線、全非零K、全higher sourceで c因子又はPf零を得る。

## 7. 範囲と次の残余

三つの全二parameter normalSS families は、旧G次元上界からではなく
actual shearと全source compatibility／factorで pointwise 排除された。
対象のm二では全P0に独立T,Vがあり、剛性後のP0Pfも非零の非被約因子型である。
全first-symbol normalizerとfinite shearの逆で、元の全sourceへ戻る。

三parameter normalSS MIXは別に残り、本稿のgeneric shear-rigidityを適用できる全strata以外のcurve／pointをまだ認識していない。
nilB0二／三parameter、Cartan既被覆、low21、normalnil b非零、co-Higgs k一等も担当分離した別義務。
新component増分は零。元reports／code／主稿／global／Gitを変更せず、同名manifestへexact inputsを保存する。
