# cubic cone 残一軌道：大域上界の独立監査と、残る T39 分岐

2026-11-08、PRIMARY E、別作者 `exceptional_coverage`。
新しい operator、point、minor、Jacobi、jet 計算は行わない。
本稿の判定は **非正規・非plane・rank二の大域上界は適用可、無条件の成分性は未証明** である。
T39 の次元を28とする推論は使えない。公開済み `57f8e29` の台帳でも T39 は25次元である。

## 1. 元の点と、閉条件から一般点へ渡せる範囲

対象は `CUBIC_CONE_RANK_TWO_ONE_ORBIT_FRONTIER.md` の canonical 点と、その24次元 actual 軌道である。
既存の全係数証明書から、実五次式 F は絶対既約、Sing(F) の被約支持は twisted cubic cone S のみで、
singular plane は無い。projective tensor の零集合は曲線とその頂点を合わせても一次元である。
全 D の rank100、Euler24、actual M49 から projective 接次元は25、軌道次元は24である。
一本法線の二次余核が零だったことを、成分性又は smoothness の証明へ読み替えない。

rank四かつ実 Pfaffian 非零の開集合では、実五次式の projective 係数は正則写像である。
被約既約でない五次式は、有限個の正次数分割の projective 掛け算像の合併に入り閉である。
したがってこの点を含む各既約 Poisson 成分の一般点でも実五次式は被約既約となる。
Sing(F) に plane が含まれる条件は Grassmannian 上の閉 incidence の proper 像なので閉である。
projective tensor の零集合が二次元以上である条件も、projective universal zero family の上半連続性から閉である。
よって同じ一般点では singular plane も rank零曲面も無い。

被約既約非正規 hypersurface は S2 なので、R1 の破れから少なくとも一つの integral 導手曲面を持つ。
Hamiltonian 導分は正規化・導手を保つので、その一般 restriction rank は零又は二である。
上記二つの閉条件は degree一と rank零を排除する。
全 length 式 `sum(delta_i*deg(S_i))=6-g` により、残る degree は二から六。
退化した degree三以上の S が hyperplane に入る場合は、ideal `(h,q_e)` の symbolic square と degree五から
実五次式に h 因子が生じるため空である。degree二の smooth quadric と quadric cone は別に保持する。

## 2. 非plane・rank二の全像上界は、一般 beta に限った上界ではない

以下は固定 surface の上界だけではなく、surface、projection、mark、frame、moduli と真の fiber を戻した
**全 PGL5 飽和射影像**の上界である。primitive beta の像が proper という別の結論を使っていない。

| 導手曲面の degree / 型 | 全源の結論 | 今回確認した最終入力 |
|---|---|---|
| 2、smooth quadric | 実 hyperplane 因子、従って空 | `SMOOTH_QUADRIC_POISSON_CONORMAL_HYPERPLANE_FACTOR.md` とその別作者 review |
| 2、quadric cone | 全像24以下 | `QUADRIC_CONE_RANK_TWO_SEPARATE_REVIEW.md` |
| 3、smooth cubic scroll | 全像24以下 | canonical-root 排除、A25 全源、boundary review、coordinate/pole erratum review |
| 3、twisted cubic cone | 全像24以下 | `CUBIC_CONE_RANK_TWO_SEPARATE_REVIEW.md` |
| 4、全 integral 曲面 | 空 | `DEGREE_FOUR_CONDUCTOR_IRREDUCIBLE_QUINTIC_ALL_EMPTY.md` |
| 5、sectional genus零の全正規化 | 全像24以下 | `RANK_TWO_QUINTIC_SURFACE_GENUS_ZERO_SEPARATE_REVIEW.md` |
| 5、genus一 elliptic 正規化 | 空 | `RANK_TWO_QUINTIC_SURFACE_GENUS_ONE_ELLIPTIC_SEPARATE_REVIEW.md` |
| 5、canonical index二正規化 | node 全像24以下、cusp 空 | `RANK_TWO_QUINTIC_CANONICAL_INDEX_TWO_SEPARATE_REVIEW.md` |
| 6、全正規化 | 全像24以下 | `RANK_TWO_SEXTIC_SURFACE_SEPARATE_REVIEW.md` |

smooth scroll の旧 integral-R 部分ノートだけでは全源を覆わない。
最終 A25 証明では R の全重複、全 binary kappa partitions、c を因子に持つ Omega、
全 full P2 fibers、第二 row の divisorial zero が保持される。
二つの primitive root 型を canonical-root 空性と A25 全源へ尽くした後、
node の真の二 scales、cusp の二 scales と加法を引いて固定 S の像六以下、
Aut(S) 六を戻して全像24以下となる。
`CUBIC_SCROLL_A25_BOUNDARY_SEPARATE_REVIEW.md` の補正通り、
reducible Omega 段階で residue の poles を L との交点だけとしない。
coordinate/pole erratum も最終 chain に含める。

degree五 genus零では全五 scroll 型と plane 導手の全 class を保持する。
残る113/023/014の projection source は13、frame24、実 Aut は13/14/16で、
014 の追加 fiber choice 一を含め全像24/23/22。
fixed-D の全 projective tensor fiber が高々一点という留数の一意性を使う。
elliptic genus一では全 zero-summand、Pic零 twist、single/non-split/split node と cusp が扱われる。
canonical index二では全非terminal/canonical境界と dependent marks を扱い、
node の全 presentation25から実 C* fiber 一を引く。cusp の短い quadric-degree 矛盾は全源を覆う。

degree六では113/122/014/023/005の全型を保持する。
023の全有限射影族上で実導手が Hom-duality と MCM base change により変動し、
nonintegral conductor は真の閉部分集合となる。projection14、frame24から実 Aut14を引いて全像24以下。
これは一つの良い射影の固定像だけの上界ではない。

degree四の全 almost-minimal 分類、degree五・六の全正規化分類を含め、
normal、smooth、ACM、clean、ordinary node、squarefree R、vertical 成分無し、
normalization-flat といった未証明の追加仮定は今回の表に無い。

## 3. 正しい大域的な二者択一

この点を含む既約成分 Z の一般実五次式が非正規なら、§1–2 の有限個の全源の像を用いて `dim Z<=24`。
元の actual 軌道は24次元なので、その軌道閉包が Z に等しい。

一般実五次式が正規なら、既受理の normal actual divisor theorem によって一般点は T39 に入る。
T39 は閉であり Z は全成分なので、この場合は `Z=T39`。
T39 の projective 次元は **25**、S61 次元は16である。
公開 Git `57f8e29:research/COMPONENTS.json`、blob
`47a41fca488c5658c7ecaca71f9a693f659fde68` の T39 entry にもこの値が保存されている。
接次元25はこの分岐と矛盾しない。

従って現在正しい結論は、各 Z が「24次元軌道閉包」又は「25次元 T39」のいずれか、という二者択一だけである。
T39 非所属を証明していない以上、新成分又は成分総数の増分は主張しない。
別途進む nodal elliptic quintic の actual T39 認識が陽性なら、必要な E 帰属はそのまま閉じる。
T39 へ入ることと、同じ点を通る他成分の不存在は区別する。

## 4. 監査の範囲

本稿は既受理の全 source proof を今回の大域推論へ適用できるかを別作者として読解した scope 監査である。
既存の各有限行列又は分類定理を新しい独立計算として再実行したものではない。
新しい数値実行は無く、純証明の receipt と入力 hash manifest を別に保存する。
global 台帳、component count、main manuscript、Git は変更しない。
