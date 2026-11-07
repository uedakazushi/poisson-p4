# smooth CI(2,2)上で零となる全Poisson bivectorはrank二以下

2026-11-07、主義務G。S⊂P4を滑らかな二quadricのcomplete intersectionとする。
S上で零となる全global Poisson bivectorは一般rank二以下である。
従ってSはrank四Poisson構造のrank零曲面にならない。
normalization、conductor、ordinary node、vertical成分の仮定は一切不要である。
非零rank二bivectorは存在するので、全Poisson sectionが零という強い命題は主張しない。

## 1. 全sourceは二つのconstant skew matrixで尽くされる

`E=Λ²T_P4` とする。Koszul列をEでtensorすると

`0→E(−4)→E(−2)²→E⊗I_S→0`。

外積Euler列とP4のline bundle cohomologyから

`H⁰(E(−4))=H¹(E(−4))=0`,
`H⁰(E(−2))=Λ²C⁵`、dimension十。

例えば外積Euler列は `0→T(−k)→Λ²C⁵⊗O(2−k)→E(−k)→0`。
k=2,4について通常のEuler列とline bundleの中間cohomology消滅を使えば上記が従う。
従って全sectionは**一意に**

`Π=q₁A+q₂B`                                                (1)

と表され、A,Bは全十係数のconstant skew 5×5 matrixである。
これは一切片でなく全20係数の完全sourceであり、Euler商の余分なkernelも無い。
Poisson条件は全350projective Jacobi `Eul∧[Π,Π]=0` を使う。
raw strict Jacobiの零を余分に課していない。

## 2. 全smooth pencilをpairwise-distinct eigenvalueで覆う

smooth CIのpencilには非退化なquadric q₁が存在する。
実際、pencilが全て退化と仮定し最大rankのmemberを取る。
rank≤3ならそのradicalはdimension二以上であり、もう一方のquadricに関する非零isotropic vectorがある。
この点はbase locus上で一方のgradientが零となるため特異点である。
最大rank四ならradical vector vについてdet(q₁+tq₂)の一次係数がq₂(v)に比例する。
detが恒等零ならq₂(v)=0で、やはりbase locusが特異となる。

q₁を非退化に取り、q₂をq₁-selfadjoint operator Mで表す。
重複eigenvalue λがあれば、eigenspaceのdimension二以上の場合にはq₁-isotropic vectorがある。
dimension一でもJordan chainがあれば、selfadjoint性からそのeigenvectorはisotropicである。
いずれもq₁=q₂=0でgradientが従属となり、Sのsmoothnessに矛盾する。
従って全eigenvalueは相異なり、orthogonal eigenbasisと非零coordinate scaleにより

`q₁=Σxᵢ²`, `q₂=Σλᵢxᵢ²`, `λᵢ≠λⱼ (i≠j)`。               (2)

逆にこの条件ならbase locusにgradient従属点は存在しない。
従って(2)は全smooth pencilを覆う。pencilのGL₂変更は(1)の全A/B源を保つ。
λを固定した一つの有理標本から全smooth Sを外挿していない。

## 3. 全Jacobiのparity blockと三Pfaffian

全20parameterの二次単項式は210個。
coordinate sign flipに対してAᵢⱼ,Bᵢⱼのcharacterはeᵢ+eⱼである。
characterが四要素集合Iとなるblockは12列であり、次の三polynomialを全て含む。

`Pf(A_I)`, `Pf(A_I+B_I)−Pf(A_I)−Pf(B_I)`, `Pf(B_I)`。      (3)

全350Jacobi係数を元の(1)から生成し、その各rowのparameter supportが単一parity blockに入ることを保持した。
まずλ=(0,1,2,3,4)のexact QQ診断では、全五blockがrank六で(3)がrowspaceに入った。
この診断自体は全λの証明として使わない。

次にI=(0,1,2,3)について、選んだ六つの元Jacobi rowをQQ(λ₀,…,λ₄)で扱った。
その六次minorは正確に

`−64(λ₀−λ₁)³(λ₀−λ₂)(λ₀−λ₃)(λ₀−λ₄)`。

三つのtargetは六つの元Jacobiの線形結合として恒等的に表される。
全functionalsを保存し、cancel後の全分母は2及びpairwise differenceの積だけである。
従って全λᵢ相異なる点を覆い、別のminor零の補集合は無い。
S5のcoordinate permutationを各blockへ実際に運び、
multivector添字のorientation signも含めた全五blockのpolynomial identityを独立に検証した。

従って全projective Poisson条件上で、各Iについて三つの(3)は全て零となる。
全A/B parameterを保持し、A又はBのrank、非零minor、chosen leafを追加仮定していない。

## 4. 全rank四排除

任意のscalar t,sについて

`Pf(tA_I+sB_I)=t²Pf(A_I)+ts PfMixed(A_I,B_I)+s²Pf(B_I)=0`。

全五Iで成立するため、skew matrix tA+sBはrank二以下である。
特に各pointで `Π(x)=q₁(x)A+q₂(x)B` はrank二以下、projective πもrank二以下。
正しい実anticanonical quintic `Eul∧Π∧Π` も恒等零となる。

Aを非零decomposable constant skew matrix、B=0とすればq₁Aは非零rank二Poissonとなる。
従ってこの証明はrank二の既存支持を誤って除いていない。
Sをrank零曲面に持つrank四の全枝が空という、必要な弱い結論だけを使う。

## 5. 独立確認と保存物

producer `code/classification_quadric_complete_intersection_zero.py` は全350×210係数を構成し、
12列blockの全functionals、symbolic minorと全分母を保存した。
QQ診断0.31秒、symbolic証明6.248秒。

独立 `code/classification_quadric_complete_intersection_zero_review.py` はproducer/CASをimportせず、
q₁A+q₂BそのものをFraction疎多項式で微分して全350式を再構成した。
全QQ診断、symbolic六次minor、全Vandermonde分母を確認し、
全五blockの三targetを各λの多項式として分母を払って検証した。
PASS、0.191秒、最大RSS16,384KiB。

証明書は
`CLASSIFICATION_SMOOTH_CI22_ZERO_PARITY_{DIAGNOSTIC,SYMBOLIC,REVIEW}_20261107.json`。
各実行前gateは `jobs/smooth_ci22_zero_parity_{diagnostic,symbolic,review}.json`。
symbolic producerと独立reviewerの完全bytesを各証明書へ保存した。

上記Koszul/Eulerの完全sourceとsmooth pencilの全normal formは別の幾何学的証明であり、
計算PASSだけがその仮定を保証したとは扱わない。
特異CI pencil、S上で一般rank二のrestriction、他のsurface、既約非正規五次全体はこの定理に含めない。
