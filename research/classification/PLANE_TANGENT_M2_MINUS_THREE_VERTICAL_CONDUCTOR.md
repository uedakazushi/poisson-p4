# t=−3 canonical tensor：X122上でvertical、residueはdηの定数倍

2026-11-08、`/root/quartic_review`。純計算による補足、新しい演算job無し。
全source／normalizationは `PLANE_TANGENT_M2_MINUS_THREE_WHOLE_ORBIT.md` と
`PLANE_TANGENT_M2_MINUS_THREE_X122_NORMALIZATION.md` に固定する。
既知成分への平滑化／全relative descentは別担当の未受理証明に委ね、本稿だけでT40等を採録しない。

## 1. projective chartの全括弧を用いる

u=1 chartの座標(a,b,c,v)では、raw homogeneous bracketsをそのまま制限せず、
`{xi,xj}_aff=Pij−xj Pi,u+xi Pj,u`を用いる。
canonical P(0,1)から

    {a,b}=3ac+3bcv+av²,
    {a,c}=av+3bc+3c²v, {a,v}=−3(a−cv²),
    {b,c}=a−cv², {b,v}=−3(b+cv+v³/3), {c,v}=0.

A=a−cv²、B=b+cv+v³/3とすると、Leibniz展開で正確に

    {A,B}=0, {A,c}=Av+3cB, {A,v}=−3A,
    {B,c}=A, {B,v}=−3B, {c,v}=0.                      (1)

を得る。全affine Fは3(A²−vAB−3cB²)。
{A,B}=0の確認を落として、c,vの二括弧だけからvertical性を結論してはいけない。

## 2. normalization base ηがCasimir

D上のdense B!=0でη=A/Bとすると

    {η,B}=0, {η,v}=0,
    {η,c}=ηv+3c−η²=−F/(3B²)=0 on D.

η,B,vはnormalizationのこのchartの独立座標で、c=(η²−ηv)/3。
従ってinduced rank二Poissonは正確に

    π_D=−3B ∂B∧∂v,

すなわちbaseη=s/tにverticalである。
normalization mapの各P2 fiberはprojective Poisson planeになる。
raw P-normal rowがEuler tangent rowの倍である場合を非不変と誤認しない。
例えばη=0 fiber(a=c=0)はraw P_a,j=−3b E_jで、affine quotient後は確かにtangent。
一方、先行候補(c=u=0)はfull Euler-row比較でもPoisson planeではない。

## 3. horizontal conductorのpoleは明示的に相殺する

volume da∧db∧dc∧dvはdA∧dB∧dc∧dvに等しい。
F_c=−9B²なので、hypersurfaceのdualizing three-formは非零符号を除き

    Ω_D=−dA∧dB∧dv/(9B²).

normalizationではA=ηBだから

    ν*Ω_D=−dη∧dB∧dv/(9B).

上のπ_Dとcontractすると、residue one-formは非零定数倍のdη（この規約ではdη/3）。
horizontal conductor B=0でのpoleはB factorにより正確に消える。
η=s/t上のdηはnormalization P1で2(t=0)のpoleを持ち、
CP=t²のcuspidal genus-one baseのRosenlicht dualizing generatorに一致する。
全polesを無根拠にlogarithmicと仮定した結論ではない。

このexplicit中央同一性は、partial normalization Z0=P_{B0}(E0)上の
vertical conductor bivectorを平滑化する別のwhole幾何証明が得られた場合に、
special tensorが本canonical点に一致することを保証する入力となる。
ambient liftの一意性はH0(Λ²TP4(−5))=H0(Ω²P4)=0で処理出来る。
relative smoothing／finite conductorのbase-changeを本稿から推測しない。
