# m2 constant-root・nilpotent diagonal接続の全空性

2026-11-07。主義務G/E。全sourceの数学的証明案。計算jobなし。
実D被約既約quintic、rank二restriction plane、reduced nonzero R、位数二、normal quadratic generic squarefreeとする。
root degree flag(0,q)の全primitive源を対象とする。
Rがtangent conic+line又はconcurrent triangleである全constant-root sourceの空性を示す。
本稿の新しいcurl identityは別作者の独立照合を要する。

## 1. 全weight零源とnilpotent diagonal fields

constant rootをuとし、G2=κ0 h u(Aq u+Bq v)、hBq=f、κ0≠0。
root-line proofによりβはlower triangular、diagonal entriesはclosed logarithmic forms。
canonical tracefree mixed gaugeでは

\[
 P_0^{ij}=(P_f)^{ij},\quad
 P_0^{i,u}=uX_i,\quad P_0^{i,v}=uZ_i+vY_i,\quad
 P_0^{u,v}=k_{20}u^2+\kappa uv.
\]

i,jはbase三座標、X,Y,Zは全linear base fields。
βuv=0からu-columnはu倍、targetG2のv² coefficient零からK02=0となる。
off-diagonal Zは全て保持し、log性を仮定しない。

tangent conic+line f=z(xz−y²)では、closed diagonal logarithmic formsに対応するX,Yは
2y∂x+z∂yのconstant倍である。
concurrent triangle f=xy(x−y)ではX,Yは(a x+b y)∂zという全二次元spaceに入る。
両場合ともtrace X=trace Y=0、X(f)=Y(f)=0、X∧Y=0であり、全linear combinationがnilpotent matrix。

actualD2のuv coefficientは3fκ。
他に現れるX∧Y、X(f)、Y(f)の寄与は上の恒等式で全て零、Zはu² coefficientだけに寄与する。
target G2のuv coefficientはκ0f≠0なのでκ≠0である。

## 2. fullweight一の一つのbad coefficientだけを使う

完全normal filtration P=P0+P1+P2を保つ。
P1のu-columnのuで割れない項をv²C、Cはconstant base vectorとする。
他の全P1、P2は固定しない。
Jac(P)のprojective四vector E∧Jac(P)=0を使う。

strictJacのnormal v³係数は

\[
 [\operatorname{Jac}(P)]_{iuv,v^3}=\kappa C_i+(D Y\cdot C)_i,
\]

strictJacのnormal v²係数は

\[
 [\operatorname{Jac}(P)]_{iju,v^2}
 =[(2Y-X)\wedge C]_{ij}-(P_{C(f)})_{ij}.
\]

Z,k20はuを含む別normal monomialにしか寄与せず、ここへ入らない。
higher normalweightのJac(P1,P1)、Jac(P0,P2)等もこのweight一係数へ寄与しない。
従ってprojectivequad(i,j,u,v)のv³係数は正確に

\[
 E_S\wedge(\kappa C+D Y\cdot C)
 -(2Y-X)\wedge C+P_{C(f)}=0. \tag{1}
\]

これは全sourceへの必要恒等式であり、chosen sliceではない。

## 3. curlによるinvertible scalar+nilpotent

base三bivectorsをvolumeでvectorへ双対化し、curlを取る。
P_{C(f)}はgradient vectorなのでcurl零。
constant vectorDに対しcurl(E_S×D)=−2D、tracefree linear fieldVに対しcurl(V×C)=DV·C。
式(1)から

\[
 (2\kappa I+4D Y-D X)C=0.
\]

4DY−DXは上の全nilpotent spaceに入り、κ≠0なのでこのmatrixはinvertible。
ゆえにC=0。

## 4. 全actualDの不可約性との矛盾

canonical全30/9源でC=0なら、P1のu-columnもu倍であり、P2にはnormal columnsが無い。
normal bracketはu(k20u+κv)なので、fullPの全u-rowがuで割れる。
u=0でprojective tensorはhyperplaneの三次元tangent spaceに入りrank≤2、
従ってactualDはuで割れる。実D既約に反する。

これでtangent conic+line及びconcurrent triangleの全constant-root m2squarefree sourceが空となる。
nilpotent fieldsが零になるperiod/resonanceも同じinvertibilityに含む。
κ0=0はD2零でm三以上の別義務であり、無断で補集合を落としていない。
新しい350×30 profile、minor、GB又はjetを一切使わず、一つのglobalcurl恒等式だけで閉じる。

triangle/transverseconicの全constant-root源は別noteでactualSNC→MPS。
flag11は別conic/split proof、flags03/12はline-class矛盾。
これらの独立receiptと本curl proofの独立確認が揃えば、reducedRの全m2squarefree sourceがMPS又は空となる。
normal double quadratic、nonredR、m3等は未完了。
