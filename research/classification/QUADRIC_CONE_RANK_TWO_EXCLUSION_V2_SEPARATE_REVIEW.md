# Rank二quadric cone全排除V2：全Fu入口と全極条件の独立監査

2026-11-08。PRIMARY E。**V2全源純証明PASS**。
対象は `QUADRIC_CONE_RANK_TWO_FULL_BASE_REDUCTION_V2.md`（9159f0fd…）と
`QUADRIC_CONE_RANK_TWO_IRREDUCIBLE_QUINTIC_EXCLUSION_V2.md`（5948af0e…）。
新CAS、point、operator、kernel、rank、minor、profile又はjetは実行していない。

結論は、被約既約actual quinticを持つrank四Poisson tensorの非正規支持に
一般restriction rank二のintegral quadric coneが現れる全枝は空、である。
可約Dの全component帰属又は全研究の完了は主張しない。

## 1. V1の欠落を復活させた全sourceだけを入口とする

旧独立全33／Euler商／actual normalizer／Hamiltonian pencil証明を入力とする。
Dpencil零は元Jacobiの帰結として使い、初めからOreと仮定した旧失敗sourceは使わない。
`ENTRANCE_ERRATUM`の全Bを元PENCILの全34 direction配列と照合し、
SSのp26,p27,p30,p32及びnilのp26,p27,p29,p30を含むことを確認した。
V1のsix-base全入口・whole結論は不採用のままbytesを保持する。
旧pole／weight計算を使うのは、V2の実正常化がその式の仮定を得た後のconditional lemmaとしてだけである。
C零は別に独立受理した全source factor証明（86924013…）をそのまま使う。

全源はP=B+(Y+vT)∧∂v、degY≤3、Y3=−qE、q=ac−b²。
Faff=G+vK、G=bvec·Y、K=bvec·T、column divergenceをZ_Bとする。
規約 `P gradF=−F divP` から

    B#dK=−K Z_B,
    B#dG+K Y=−G Z_B+G T

を独立に導いた。K非零ではY=gT−B#dg、g=G/Kである。
K零へこの除算を流用せず、各K零枝の全polynomial operatorを別に使う。

## 2. SSの全Fu-only、分母零及び追加二枝

全JacBのdegree二を再展開すると、V2の
s零、(κ+1)ξ=(κ+1)ν零、κt+(1+2κ)ζ零、κr+η零を得る。
κ非−1かつw=1+2κ非零で、実a′=a+r,c′=c−t/wにより
canonical Bへ戻ることを手で確認した。d′=d−trである。
w零（κ=−1/2）はt零、η=r/2、ζ自由を保持し、a′=a+r,c′=c−2ζで処理する。
各変換はglobal affine/PGLで逆は負translation、全Yもpushする。
残JacBがA=D零、Lie_TBが全T0零とρd零を与える。κ=−1をこの枝へ捨てていない。

κ=−1では全B=Jac(f)にquadratic terms2ξa²−2νc²が残る。
Lie_TBの全degree二からξ(2−σ)=ν(2+σ)=0及び
T0=((1+σ)r,0,(σ−1)t)を得る。
a′=a+r,c′=c+tで全T0とab/bc項が消え、下位constantを新A,d,Dへ可逆に吸収する。
degree一の条件も全σで保持する。
したがって旧canonical cubic枝に加え、必ず

    σ=2: f=λa²+b³−abc, T=(3a,2b,c), λ非零;
    σ=−2: f=λc²+b³−abc, T=(−a,−2b,−3c), λ非零

の両方がある。λ零は旧canonical枝へ戻る。

第一枝はK=6f、primitive c-linear fは絶対既約でfield C(a,b)。
f因子が無ければη=G/6 mod fは非零polynomial classで、全pole cancellationから
B#dη=ηT。B_ab=−abによりaη_a=2η,bη_b=−3ηなのでη=C a²/b³。
c一,b=t,a=t(1+√(1−4λt))/(2λ)はf上の正しいalgebraic formal branchで、
a,bともvaluation一、従ってηのvaluationは−一。
polynomial classのpullbackがC[[t]]に入ることに反する。
これはPoisson jetの近似でなく、一つの正確なDVR上の代数関係である。
第二枝もa一、cを同じquadratic rootとしてη=C c²/b³のvaluation−一で排除する。

## 3. nilの全centerとell/hの保持

全JacBのdegree二から
t=(12−8κ)ξ、(2κ+3)(κ−1)ξ零、κr=(4κ−2)ζ、κs=(3−4κ)ηを再導出した。
m=1−κ、R=r−2ζ、S=s+η、L=t−4ξなら
Z_B=2(R−2mb,S−mc,L)。

m非零／ξ零は、κ零を割らずb′=b−R/(2m),c′=c−S/mでcenterする。
divergenceの全constant零とJacB、すなわちbvec·N零の全係数比較から

    B=B2(κ)+ell Jac(q)+h c∂a∧∂b+A∂a∧∂b,
    T=σE+wN+τ∂a, w=1−2κ,
    σell=0, 2κτ−σh=0, ellτ−2σA=0

が**全**必要形となることを確認した。N=2b∂a+c∂b。
ell/hを旧six表示に勝手に落としていない。
σ非零ではell=A零、κ非零はa′=a+h/(2κ)でh,τを同時に消す。
κ零ではh零、a′=a+τ/σでτを消す。
σ零／κ非零ではτ零、同じa-shiftでhを消し、ellと変化後Aを全て残す。

σ零で残る全K零源のoperatorは

    [(c+ell)∂c+(wq−A−h c)∂q]G=(1+2w)G,
    N(G)=0、従ってG∈C[c,q]。

κ非零ではhを除いたU=c+ell,Q=q−A/wのweight比較を使う。
w零も別のpolynomial比較でU因子となる。
U因子を避けるκ=3/4はG=α−(5/2)U(q+2A)²だけで、
bvecがproper ideal(U,q+2A)に入るためα零となる。
quotient C[a,b]/(−ell a−b²+2A)はell=A零でも非零なので、隠れたunitは無い。

σ=κ=τ零ではU=c+ell,V=q−A+h ellとするとoperatorは
E_(U,V)−hU∂V、固有値三である。
total degreeごとに定数＋nilpotentなのでdegree三だけが残り、h非零ならU³だけ。
h零でもspatial degree≤5がV³を禁じ、全GにU因子がある。
τ非零ではell零、K=τc²。全hを保つpole係数η=G|c0/τはdegree≤5のpolynomial。
全tangent derivativeを含む二成分からη_a零、(b²+A)η′=(6b+τ)ηを得る。
非零polynomialの最高degreeは六を強制されるので矛盾する。

## 4. nil κ一と新κ=−3/2・ξ非零を落とさない

κ一は全quadratic Hamiltonianを保持する。
σ非零ではT1=σI−Nが可逆で全T0を実translationにより消し、
degree二／一上の−σ−N、−2σ−Nが可逆なのでf2=f1零となる。
σ零では全係数比較からξ零、T0=(−2ζ,3η,0)。
実b′=b+ζ,c′=c−3ηとa-shift後、全f=(c+ell)q+A c、T=−Nを残す。
weight式p−j=−1、spatial degree≤5から全Gにq+A因子がある。

新κ=−3/2／ξ非零では、[T,Z_B]のc成分がσ零を強制する。
残りからT=−4Z_B/5+δ∂aで、Lie_ZB B零と∂aBabのc係数−三からδ零。
JacBはbvec·Z_B零なのでK零。
modular恒等式にdGを掛けると、G非零ならZ_B(G)零。
center後のZ_B=(-10b,−5c,40ξ)と、固定最高部G5=2c(ac−b²)²を独立に確認した。

## 5. Jordan4の全degree≤5 invariantを手で尽くす

Z_Bをhomogenizeし、非零diagonal rescalingだけで
D=h∂x1+x1∂x2+x2∂x3とする。x1はcの非零倍であり、最高部のa²c³は消えない。
degree五空間Sym⁵(Sym³ C²)のweight multiplicitiesをpartitionで数えると、
最高weight15から下のk=0,…,7について1,1,2,3,4,5,6,6。
差分は最高weight15,11,9,7,5,3を一つずつ与え、dimensions16+12+10+8+6+4=56で全体を尽くす。
char零のsl2完全可約性によりkerDは六次元である。

I=2hx2−x1²、J=3h²x3−3hx1x2+x1³とするとD(I)=D(J)=0。
次の多項式を直接展開した。

    Δ=(J²+I³)/h²
      =9h²x3²−18hx1x2x3+6x1³x3−3x1²x2²+8hx2³。

h⁵,h³I,h²J,hI²,IJ,hΔは全てdegree五、D零で、上の六つのweightが相異なるので完全基底。
h零へ制限すると残るのはIJ=−x1⁵だけである。
従って任意のaffine degree≤5 invariantの最高degree五部はc⁵のconstant倍に限られる。
固定G5のa²c³係数二とは両立せず、新whole枝も空となる。
個々の不変式例だけから完全性を推論したものではない。

## 6. 残る旧conditional pole lemmaとglobal因子

SS canonicalのb-pole residue及びnil canonicalのc-pole residueは、分母とbvecの消去後に
Y_b／Y_cのpolynomialになる。分子のtangent derivativeも含めてpole式を使う。
SS κ=−3/2のq-poleではη=(∇q·Y)/(4ρ) mod qがgraded polynomial classであり、
全pole式がEuler(η)=−4ηを強制するため零となる。
SSの二unimodular cubicではη²(1−ab)又はη²(1−bc)が非零constantとなるが、
function fieldの対応prime valuationが奇数で不可能。
その因子がaffine surface上unitでも、field valuationの議論は有効である。
これら旧代数lemmaはV2の全入口正常化後に限り適用した。

各proper chart factorはdegree一／二／三で、cofactorのdegreeは五から差し引いた上限を持つ。
homogenizationで元actual degree五にもproper factorがあり、cofactor零ならrank四でない。
Euler gaugeと実normalizer逆はこの結論を変えない。
C零／SS／nil、全κ,σ、全Fu-only、全unit零補集合及び全Y,T0を合わせ、冒頭の全空性を受理する。
V1の誤入口を事後的に正しかったと呼ばず、明示withdrawalとV2新identityを保存する。
