# 固定319四次行列の欠損は余次元二以上

2026-11-03。これは新しいconceptual transferである。旧H8 transverse theoremが証明した「何らかの適格cutで長さ65」と、固定weighted chartの全四次行列との間を、nilpotent thickening上のparameter changeで結ぶ。alternate cutと固定cutのrow spacesが等しいとは仮定しない。新しい係数行列・GB・field・点計算は行わない。

## 1. 正確な範囲と入力

baseは U0=Spec QQ[a,e,d,1/(aFBC)]。nuisance=0の同じH8表示に限定する。対象M4は、HYPERPLANE_WEIGHTED_CHARTの固定八normal variablesによる、元135式のweight4 ideal componentである。354×239表示と既存319×239表示は同じimageを持つ。T34または既知unit diagonalで移したTnを**全行で一貫して**使う。最近のoracle/frame adapterを省いて混在させた行列は対象でない。

依存する既受理結果は次の通り。

- TRANSVERSE_CODIMENSIONおよび独立REVIEW: A³coreにQQ開U_anyがあり、補集合の次元≤1。各点で有限個の適格cutの一つによる横断上界65があり、その同じH8 relative completion Sは、実在binary部分の正則completion R（H8 base8+binary3=11変数）上で自由65。restriction R→S→Rが分裂し、同じbinary componentのgeneric local length65が全rankを使い切る。
- WEIGHTED_CHART: 固定cutは実在binary三方向への射影determinant−C/25を持つ。34選択式、正の八normal weights(1,1,1,1,1,2,1,1)、exact degree≤4 implicit graphが全U0上で定義される。
- QUADRATIC_CODIMENSION: weight2 rank<6のcore集合Z2は次元≤1。
- CUBIC_CODIMENSIONと既存全点upper57: weight3 rank<57のcore集合Z3は次元≤1。

U_anyは固定cutそのものの全射開として定義されていない。特に旧C-special stepはalternate cutを使う。ここを以下のlemmaで補う。

## 2. Nilpotent thickening上のcut変更lemma

Bをcomplete regular local ring、R=B[[s1,s2,s3]]とする。local complete algebra SにR→S→Rのretractionがあり、SはR上自由n、restriction primeにおけるgeneric local lengthもnとする。

N=ker(S→R)と置く。K=Frac(R)上でdim S_K=n。このdimをrestriction primeのlocal length nが全て使うため、S_Kは一つのArtin local factorだけを持ち、N_Kはnilpotentで(N_K)^n=0。SはR-torsion-freeなのでN^n=0。S/N=Rはreducedだから、S_red=Rである。

f1,f2,f3∈m_Sについて、そのreductionsがB上のRのregular coordinatesであるとする。centered cutsを使い、R'=B[[z1,z2,z3]]→Sをzj↦fjで定める。Rへの合成はformal inverse function theoremにより同型である。

N^j/N^(j+1)はfinite R-moduleで、R' actionはS/N=Rを通る。R'≅Rなので各quotientはfinite R'-module、N^n=0のfinite filtrationからS自体もfinite R'-moduleとなる。

Sは旧R上freeなのでCohen–Macaulayである。新R'は同じdimensionのregular local ringで、Sはfinite local algebra。従ってSはR'上maximal Cohen–Macaulay、Auslander–Buchsbaum（またはregular system of parameters）によりfinite freeとなる。そのrankは唯一のgeneric componentに沿うintrinsic length nであり、R'→S_red=Rが同型だからnである。

したがって **どの適格cut fも同じ長さnを与える**。S/(f1,f2,f3)はB上自由n、その閉fiber長もn。これはcutのrow spacesの等号ではなく、同じnilpotent厚みのbase parameter changeである。Sが旧R上free、restriction componentが全generic rankを使うこと、cutsのreduced-coordinate条件のいずれも省けない。

## 3. 現在の固定cutにlemmaを適用する

固定cutは、生の古いindex(5,19,23)との同一視で定義しない。WEIGHTED_CHARTのchart_basis Wを使い、W^-1*yのfree34座標のうちfamily_parameters=[1,2,5]に対応する三つをfとする。既存producerはW^-1*T34のfree34 blockがI11、W^-1*U34の同blockが零であることを使っている。従ってf=0は正確にT34のnormal34八列によるfixed normal graphである。

実在binary三方向のこのfへのJacobian determinantは、同じproducerが確認した−C/25。U0ではunitなので、lemmaのnew-cut仮定を満たす。cutのsection valueは必要ならH8 base functionを引いてcenterする。closedcore fiberの固定cutは変わらない。

U_any∩U0の任意の点では、旧定理から得られるSをn=65でlemmaへ適用する。これはH8 relative completionの話であり、nuisance零のbeta自体がgeneralとは言わない。H11や全marked stratumへ拡張しない。結論は、選択されたcutがalternateでも、今回の固定eight-variable cutのcompletionも長さ65を持つ、ということである。

C=0はU0から除外されている。一方、未知height1 factorがmodpでC点へ衝突する旧議論は、U_anyの存在を証明するために既に処理済みである。今回新しくCをunitとして全A³へ戻したり、C-special pointで固定M4を全射と主張したりしない。

## 4. 低weighted Hilbert部分と全四次

点がU_any∩U0にあり、Z2,Z3の外とする。exact normal graphのpositive graded algebra Aのweight0,1,2,3部分のdimは正確に

`1, 7, 29−6=23, 91−57=34`。

合計は65である。selected graph後の全残余式はweight2..6の正確なhomogeneous polynomialsであり、完成原点局所環はこのAのcompletionに一致する。

positive graded algebraの非零homogeneous elementはgraded maximal idealでlocalizeしても零にならない。annihilatorがhomogeneousで、非零constantを持つannihilatorなら全体が零になるためである。さらにNoetherian local ringはそのcompletionへinjective。従って上の65次元部分に加えてA4の非零elementがあればcompletionの長さは66以上となる。固定cutの長さ65に反する。よってA4=0。

weight5/6 generatorsはpositive gradingのweight4 componentへ寄与できない。従ってA4=0は、元のfull weight4 matrix M4の239方向の全射と正確に同値である。G6を省かない。旧・新の双方がspecialize(...,extra=True)による同じ135式を使い、weight3のG6行とその7prolongationsも保持する。homogenized関係t*G6h=ΣK(q)F4hをunit G6除去に取り違えない。

## 5. 固定欠損のcodimension

以上からset-theoretically

`{rank M4<239} ⊂ (A³\U_any) ∪ Z2 ∪ Z3`

がU0上で成り立つ。右辺の各集合は次元≤1なので、固定M4の欠損は次元≤1、すなわち**余次元≥2**である。これはfull cubic57だけからの四次推論ではない。旧relative free65と適格cut変更lemmaが追加の本質的入力である。

全点でのM4全射は結論しない。codim2の曲線・点が残り得る。既存122 circuitのDeltaH/DeltaQ chart上では、同じrank identityによりC122 rank<64にもheight1 componentがない。

## 6. 80 directionsに対する意味と未解決gate

仮に319行sourceに対する80本の恒等kernel directionsが**各点で独立**と証明できれば、そのspanはsubbundleでquotient Eはlocally free rank239となる。generic full rankからE→O^239はgeneric isomorphism。局所square determinantのnonunit零点はheight1 componentを持つため、上のcodim≥2に反する。従ってその場合は全域全射となる。

80本がabstract module basisであることだけでは十分でない。例えば最後のimage directionが(x,y)というidealなら、generic full rankとcodim2 defectは両立し、kernelはabstract freeでもsource内での特殊化独立性を失う。必要なのは80 source directionsのpointwise rank80または同等のsplit/subbundle証明である。このgateは未証明であり、fullquartic、長65の全U0成立、ambient除去はなお未解決である。
