# normal CI22の全tangent tensor：global Euler liftの補題

2026-11-07。`NORMAL_CI22_RANK_TWO_EXCLUSION.md` §1のglobal sourceの詳細。
S=(q1,q2)⊂P4をnormal integral非退化CI22とし、πをSにtangentな
H0(Λ²TP4)のsectionとする。結論はrankを仮定しない全source lemmaである。

## restrictionの全五係数

adjunctionでωS=OS(−1)。normal性によりSregの補集合はcodimension二であり、
π|SregはωS^-1|Sregのsectionとして一意にH0(OS(1))へ延びる。
CI22のKoszul列からSはACMかつlinearly normalなのでh0(OS(1))=5。

定数五vector volumeを固定し、各linear functional ellに対して
U=ι(ell)volume、JU=ι(dq2)ι(dq1)Uと置く。
exterior contraction identityから

`E∧JU = ell(E) ι(dq2)ι(dq1)volume + terms divisible by q1 or q2`

となる。二つのnormal gradientsが独立なSregで最後のthree-vectorは非零、
従ってprojective tangent bivectorのrestrictionは非零constantを除きell(x)である。
全ellはH0(OS(1))を尽くす。所望のJUを一つ選ぶとπ−JUはSregで零、
ambient vector bundleのS上sectionはtorsion-freeなのでS全体で零となる。

## Euler商からglobalな全二十係数へ

homogeneous quadratic raw lift ΠとJUの差をΔとする。
restricted exterior Euler列

`0→TP4|S→Λ²(C5⊗OS(1))→Λ²TP4|S→0`

によりΔ|SはTP4|Sのsectionを表す。
さらに通常のrestricted Euler列

`0→OS→C5⊗OS(1)→TP4|S→0`

を使う。ACMからH1(OS)=0なので、そのsectionは五つのglobal OS(1)section WSへliftする。
linearly normal性により各成分はambient linear formへliftし、global linear vector field Wを得る。
従って同じ**global Euler gauge**を引いたΔ−E∧Wは、raw全十成分がS上で零となる。
各成分はquadratic、I_S(2)=span(q1,q2)なので

`Π = JU + q1A + q2B + E∧W`

であり、projective tensorとして全JU五係数＋constant skew A,B全二十係数に入る。
S上だけで定義されたgaugeをglobalにliftしたと黙って仮定していない。

JU自身は両qをstrict Casimirとし、q1A+q2Bのnormal contractionsはIに入る。
従ってこのraw sourceは両normalに同時tangentである。
局所generators qi/Xj²のbracketのprojector補正はI²に入り、
abelian conormalの必要条件をstrict{q1,q2}∈I²として計算出来る。
prime CIのgr_I=(R/I)[Q1,Q2]からI²はI-primaryであり、dense chartの条件は全homogeneous条件へ戻る。

この補題はnormal性を用いたcodim二延長を含むため、non-normal CI22にはそのまま適用しない。
