# 滑らかなquadric上rank二・零曲線integralの全枝は空

2026-11-07。主義務G。producer `root`。
別作者監査前の全source幾何証明。有限標本又は新しいCAS計算を用いない。

**命題。** P4上のrank四Poisson tensorの実退化因子Dが既約被約quinticであり、
その非正規曲面成分Sが滑らかなquadricとする。
S上のrestrictionは一般rank二で、その零divisor Rが被約既約と仮定する。
この全枝は存在しない。
normal coneのordinary-node性、vertical不存在、Dの正規化又は導手のsmoothnessは仮定しない。
Rの可約・非被約型、quadric cone又は他のSは明示補集合である。

## 1. 全normal symbolとparallel tensor

Sは一つのhyperplane内のquadricで、idealはI=(h,q)、deg h一、deg q二である。
prime complete intersectionなのでsymbolic square/cubeとordinary powerは一致する。
もしDのS沿い位数が三以上なら、degree五のI三乗の元は全てhを因子に持つ。
q三乗はdegree六であるため、h無しの項を作れないからである。
D既約に矛盾するので、位数は正確に二。

H=O_S(1)、S=P1×P1、H=(1,1)、K_S=-2H、
N=N_{S/P4}=O_S(H)⊕O_S(2H)である。
normal quadratic symbolを

`G∈H0(S,Sym²N*⊗O_S(5H))`

と書く。Gは零ではなく、全normal-order三以上の項もDに残している。
SのPoisson sectionをsigmaとすると、R=div sigma、R∈|2H|はintegralである。

`IRREDUCIBLE_NONNORMAL_RANK2_POISSON_MODULE.md` の全局所計算から、
{I,I}⊂I²とconormal flat Poisson connectionを得る。
U=S\Rではsigmaが非退化なので通常のflat connectionである。
Hamiltonian flowはpi二乗とsigmaを保存する。
associated gradedへ取り、adjunctionでdet T_P4|S=det T_S⊗det Nを用いると

`gamma=G/sigma ∈ Sym²N*⊗det N`

がUでparallelとなる。
Hamiltonian fieldsはTUを張るから、全方向のparallel性である。
従ってgammaはUの全点で非零、そのrankも一定である。
Gのdivisorial vertical zeroはRだけに支えられる。
normal connectionのlog性又はS上のnormal bundleの自明化は仮定していない。

## 2. rank一normal quadraticはPicard parityに反する

generic rankが一の場合、関数体上でGは一つのlinear formの平方である。
そのunique root lineをN*内でsaturateすると、滑らかなsurface上のline bundleからのmapを得る。
従ってあるline bundle Aと、divisorial common factorを持たないnormal form ell:N→Aについて

`G=gamma0*ell²`, `div gamma0=V≥0`, `2A+V=5H`

と書ける。codimension二のbasepointは許している。
このfactorizationは各codimension一のvaluationで原始rootを取り、line bundleのHartogsで延長すればよい。
§1よりVはRに支えられ、RがintegralなのでV=kR、k≥0は整数。
従って

`2A=(5-2k)H`。

Pic S=Z²、H=(1,1)なので、右辺の両成分は奇数。
line bundle Aが存在することに反する。
normal quadraticがsquareとなるcusp・より高い横断型をここで全て含む。

## 3. rank二normal quadraticの全vertical分岐

generic rank二ならgammaのdiscriminantはparallel scalarで、U上の非零複素定数である。
従ってGのdiscriminantは非零定数倍のsigma²。
関数体上で二rootは分裂する。両rootをsaturateし、primitive maps ell_i:N→A_iとして

`G=gamma0*ell1*ell2`, `div gamma0=V≥0`, `A1+A2+V=5H`

を得る。
delta=det(ell1,ell2)はA1+A2-det Nのsectionである。
discriminant恒等式から、非零constantを除き

`gamma0*delta=sigma`、従って `V+div delta=R`。

これはclassだけでなく実sectionの恒等式である。
Rがintegralかつreducedなので、全可能性はV零又はV=Rの二つに尽くされる。

### V零

blowupのexceptional E=P(N)上で、全normal-cone divisorは二つのhorizontal成分Gamma1,Gamma2だけ。
generic nodeにおけるinverse-formのclosed residueは、それぞれ非零constant lambda,-lambdaである。
したがってlog-residueのGysin関係は[Gamma1]=[Gamma2]を与える。
Gammaの特異点は解消して押し下げればよく、全点SNCとは置かない。
projective-bundle Picard公式からA1≡A2、他方A1+A2=5H。
Pic S=Z²の数値同値は整数座標の一致なので、2A1=5Hというparity矛盾を得る。
これは既受理 `IRREDUCIBLE_NONNORMAL_RANK2_CLEAN_SURFACE.md` のGysin部分の直接の適用である。

### V=R

deltaは全Sで消えず、N→A1⊕A2は同型となる。
Nの全line-bundle direct summandはO(H)又はO(2H)である。
実際、End Nのglobal matrixは、Hom(O(2H),O(H))=0とH0(O_S)=Cから三角行列。
rank一idempotentの二diagonalは一と零で、unipotent conjugationによりdiagonal idempotentへ直る。
そのimageとquotientは上記二summandだけである。

従って一方のA_iはH。そのnormal form N→O(H)は定数倍のdhしかなく、
Gはnormal coordinate hで割れる。
一方、prime CIのsquareとdegreeから実五次式は全て

`F=h² A3+h q B2+q² C1`

と書ける。
normal quadraticがhを因子に持つことはC1|S=0を意味する。
Sを含むlinear formはhの定数倍だけなのでC1=c h。
従って実F全体がhを因子に持ち、既約性に矛盾する。
高normal-orderの項を無視した議論ではなく、上の表示はdegree五のI²全体である。

## 4. 結論と残余

generic normal quadraticのrankは一又は二の全場合で矛盾した。
従ってsmooth quadric・rank二restriction・integral reduced anticanonical零曲線の全枝は空。
conormal parallel tensor、primitive root saturation、実vertical divisorとPicard parityだけを使い、
normalization・導手の全点normality、simple-poleなvertical residue又は未確認のflat familyは使わない。
R可約・非被約、quadric cone、他のdegreeやrank、全G/Eはこの命題から閉じない。
