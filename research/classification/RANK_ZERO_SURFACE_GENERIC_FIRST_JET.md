# plane以外の非正規曲面：multiplicity二とrank零Frobenius first jetの二型

2026-11-07、主義務G。D⊂P4を既約被約quintic、
S⊂SingDをintegralな二次元成分とする。
以下はgeneric normal symbolの制限であり、全Poisson支持又はG非優越を証明するものではない。

## 1. plane以外ではgeneric transverse multiplicityは正確に二

もしDがSの一般点でmultiplicity三以上なら、一般p,q∈Sを結ぶline上で
Dの五次方程式はp,qのそれぞれに三重以上で零となる。
degree五の一変数多項式は相異なる二点で合計六重に零とはなれないので、
このsecant line全体がDに入る。
従ってclosed secant variety Sec(S)⊂Dとなる。

SがP4をspanすれば、Sec(S)=P4であることを次のelementaryなjoin argumentで証明する。
projective varietyのvertex locusはlinearである。
実際、affine coneがdirection v,wのtranslationで不変ならv+wのtranslationでも不変であり、
vertex directionはvector subspaceをなす。
Sがnonlinearなら、一般p∈SはS自身のvertexではない。
従って `join(p,S)` は三次元である。
もしdim Sec(S)≤3なら、irreducibilityから `join(p,S)=Sec(S)` となり、
一般p∈SはSec(S)のvertexとなる。
linearなvertex locusはSのspan P4全体を含み、Sec(S)=P4となってdim≤3に矛盾する。
よってdim Sec(S)=4であり、Sec(S)⊂Dは不可能である。
このP4内の補題にはSeveriのhigher-ambient分類を入力しない。

SがP3をspanすれば、そのsecant varietyはP3である。
Sec(S)=SならS自身がlinearでなければならず、これはspan P3と両立しないからである。
従ってDがそのhyperplaneを含み、既約五次という仮定に反する。

S⊂SingDなのでmultiplicityは少なくとも二であり、
**span dimensionが三又は四ならgeneric multiplicityは正確に二**となる。
spanがP2となるplaneだけは、この証明ではmultiplicity三以上を排除しない。
実際、既知のtriple-plane例と矛盾しない。

このmultiplicity二はnormal quadraticの二相異なる根を保証しない。
generic transverse nodeとsquare tangent coneの両方を保持する。

## 2. rank零の一般点における四次元Lie代数

π|S=0とし、Sのsmoothな一般点pを取る。
Poisson tensorのlinearizationは `L=T_p*P4` 上のLie bracketを与える。
πの全係数はS上で零なので、その接方向微分も零であり、
`N*=N_p*(S/P4)` は二次元ideal、商 `L/N*=T_p*S` は二次元abelianである。

§1によりDのnormal二次symbolは非零である。
π自身はSで零なので、このsymbolはlinearized PoissonのPfaffianそのものである。
従ってLのlinear Poisson tensorはgeneric rank四、すなわちLはFrobenius Lie algebraである。

## 3. conormal idealはabelianでなければならない

二次元nonabelian Lie algebraはaff1であり、centerは零、全derivationはinnerである。
もしN*=aff1なら、商の二つのliftからinner derivationを引くことで、
両liftをN*とcommuteさせられる。
それらの相互bracketはN*に入り、N*全体とcommuteするのでcenter零から零となる。
従って `L=aff1⊕C²` となり、そのlinear Poisson rankは二以下。
Frobenius性に反する。
従ってN*はabelianである。
このgeneric消滅はbundle bracketとして被約smooth locus全体へ延長するので、
そこで `{I_S,I_S}⊂I_S²` となる。

## 4. Frobenius二型を全て尽くす

N*がabelianなので、quotient C²の作用はcommuteする二つのmatrix M0,M1∈End(N*)である。
商liftの相互bracketによるcocycleはあるが、Pfaffianには寄与しない。
normal covector nにおけるPfaffianは、符号を除いて
`det(M0^t n,M1^t n)` である。
非零であるため作用像は二次元で、non-scalar matrixを含む。

二次元空間のnon-scalar matrixのcentralizerは正確に `C*identity+C*M` である。
従って作用像はこの二次元algebra全体であり、quotient basisを取り直してM0=identityと出来る。
商liftのbracketがc∈N*なら、二番目のliftからcを引くことでcocycleも消える。
残るnon-scalar Mは、複素数体上で相異なる二固有値又は一つのJordan blockだけである。

したがって一般点の全first jetは、Lie代数同型まで次の二型だけである。

| 型 | 非零bracketのnormal form | normal Pfaffian |
|---|---|---|
| semisimple | [a,u]=u, [b,v]=v | u*v |
| dual numbers | [a,u]=u, [a,v]=v, [b,v]=u | u² |

前者はaff1⊕aff1で、normal quadraticは二相異なる根を持つ。
後者は `C[epsilon]/(epsilon²)` のaffine Lie algebraで、normal quadraticは平方である。
よってrank零・plane以外のSでは、generic ordinary doubleとsquare normal coneという
二つのfirst-jet機構に構造的に有限化される。

## 5. global化でまだ必要なこと

これは各一般点のLie algebra型の分類であり、S全域での一つの座標又は一つのsourceを与えていない。
二eigenlineのmonodromy、discriminant曲線、vertical共通因子、
Sの特異点、及びhigher normal係数は引き続き保持する。
smooth locus上では、abelian conormal idealに対するquotient作用を
`θ:N*→N*⊗T_S` というintegrable co-Higgs fieldとして記述出来る。
そのglobal spectral data又はnilpotent line subbundleを処理することが次の構造課題である。

特にこのfirst-jet分類だけからclean条件、q(S)>0、既知成分への全被覆、
又は未知G成分の非優越を推論してはいない。
