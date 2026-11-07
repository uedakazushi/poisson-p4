# rank二plane・位数四の全zero cubic/全normal-root支持はG非優越

2026-11-07、主義務G。全支持を結ぶproducer幾何証明。
旧独立受理済みpacketと新nonred squarefree全六sourceを合わせる。
新six-sourceの別作者receiptを受理してからaggregate claimを採録する。

**命題。** actual Dが被約既約quintic、S⊂SingDがrestriction rank二plane、S沿い一般位数が四である
全projective Poisson sourceの像は、次元二十四以下の有限個の閉族の和に入る。
πSの零cubic Rの被約性、normal leadingのsquarefree性、接続のlog性又は一意性を仮定しない。
従って未知β優越成分の一般元はこの全枝に属さない。
全点の成分帰属、m二/三、全rank二plane又は全Gの尽くしは主張しない。

## 1. zero cubicとnormal binary quarticの全分岐

πSがrank二なのでその零divisor Rはnonzero cubic。
複素plane cubicの全型はirreducible reduced三型（smooth/node/cusp）、
reduced reducible四型（transverse conic+line、tangent conic+line、triangle、concurrent三line）、
nonreduced二型（2L+M、3L）である。
Rがnonreducedなら重複成分はlineでなければdegree三に入らないので、この二型で尽くされる。

normal leadingはG4∈H0(O(1))⊗Sym4 C2、G4非零。
generic binary root partitionはsquarefree、211、31、22又は4の五つである。
coefficient image Wのrank一も全て保持する。
generic partitionを基準に全linear spaceを分けるので、chosen normal coordinateや一つの特別quarticから
whole supportを推論しない。

## 2. irreducible Rと全image-rank一

irreducible Rのsmooth/node/cusp全型は
`PLANE_IRREDUCIBLE_CUBIC_RANK_TWO_REVIEW.md`の全conormal/全30/9証明によって空。
普通二重点、leading squarefree又はsimple-pole log性を追加仮定していない。

全Rでimage-rank一ならG4=L*g4なのでD=L*g4+g5はexactly essential-three。
`ESSENTIAL_THREE_QUINTIC_EIGHT_COMPONENTS.md`とその全source別作者reviewから、全源は八23成分に入る。
これはroot partitionに無関係なpointwise recognitionである。

## 3. squarefree leadingのreduced/nonreduced両支持

reduced reducible Rの全squarefree sourceは既存の三packetで尽くされる。

- `PLANE_REDUCIBLE_CUBIC_M4_SCALAR_BOUND.md`及び別作者SCALAR review：三nonconcurrent型のgeneric j又はrank一、全像24/23/23以下。
- `PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_BOUND.md`とSPECIAL別作者review：I0/J0の全stable pencils、三nonconcurrent型で同じ上界。
- `PLANE_CONCURRENT_M4_SPECIAL_EXCLUSION.md`と`PLANE_CONCURRENT_SPECIAL_SEPARATE_REVIEW.md`：concurrentの全special rank二は共通partial因子で空。
  既rank一/generic jの全Saito矛盾も戻して全squarefree concurrentが空。

nonreduced Rではsimple-pole scalar化を流用しない。
`PLANE_M4_NONREDUCED_SQUAREFREE_RECOGNITION.md`は、I,J両非零ならdiscriminantL6を用いて
元の接続のpole supportより大きいP2\L=A2上でroot coverを分解し、G4=L*g4へ戻す。
R=2L+Mで接続がMに高位極を持ってもこの幾何は有効である。
I0/J0 rank二の全pencilsはused-line supportから六actualD正常形へ入り、
full150 canonical Casimir源の全五Pfaffiansが恒等零なのでrank四空。
従ってnonredRの全squarefree sourceはpointwiseに八essential-three成分へ入る。

## 4. repeated leadingの全支持

`PLANE_M4_REPEATED_LEADING_GEOMETRY.md`のP1 Bertiniにより、
generic repeated rootsはlinear space Wのfixed repeated rootsである。
generic22/4では全gcd degree四、従ってrank一で§2へ入る。
残る31はfixed triple root、211はfixed double rootである。

31のimage-rank二は`PLANE_M4_TRIPLE_ROOT_BOUND.md`で全Rに対して空又は全像二十一以下。
全actualD二正常形とfullcanonical150源を使い、nonunique conormal接続を一意と仮定しない。
全幾何と算術は`PLANE_M4_REPEATED_SEPARATE_REVIEW.md`に別作者確認されている。

211のreduced reducible R全支持は`PLANE_M4_REPEATED_211_BOUND.md`及び同別作者reviewで
image-rank二・三の全六profileへ尽くされ、全像二十三以下。
nonredRではcollision conic/used linesと全symmetric pencilから二sem fixed及びwholeJordan lineへ還元され、
`PLANE_M4_TWO_LINE_COLLISION_EXCLUSION.md`及び`PLANE_M4_TWO_LINE_CASIMIR_SEPARATE_REVIEW.md`により全rank四空。
Jordanの全tau、唯一walltau零、全vertex dependenceを保持したfull150 canonical源が共通partial_wを持つ。

## 5. 全支持表とGへの矢印

| 全generic normal-root枝 | Rのscope | 全sourceの結論 |
|---|---|---|
| image rank一、root任意 | 全R | 八essential-three成分、像≤23 |
| squarefree高rank | reduced reducible R | 空又は像≤24 |
| squarefree高rank | nonred R | rank四空 |
| 22/4 | 全R | rank一、八成分 |
| 31高rank | 全R | 空又は像≤21 |
| 211高rank | reduced reducible R | 空又は像≤23 |
| 211高rank | nonred R | rank四空 |

irreducible Rは全枝空として§2で含めた。
全normalizations、全base/normal/tensor parameters、全finite wallsと零補集合は各packetに含まれる。
imageのZariski閉包を取っても次元は増えない。
未知β優越成分はdimension25/26/27なので、この有限閉和の一つに含まれることは出来ない。
これが本命題からGへの正確な非優越の含意である。

位数五ならDはnormal二変数だけのbinary quinticとなり複素数体上で分解する。
従ってS⊂SingDという入口を戻すと、rank二planeの未処理重複度は二又は三だけとなる。
それらのreduced reducible及びnonred zero cubicは引き続き別義務。
本証明から全E membership又は全分類の完了を推論しない。
